#include <algorithm>
#include <cmath>
#include <vector>

#include "myplanetsim/dynamics/dry_hydrostatic_sources.hpp"
#include "myplanetsim/dynamics/shallow_water_driver.hpp"
#include "myplanetsim/dynamics/surface_orography.hpp"
#include "support/test.hpp"

namespace {
mps::ExperimentConfig williamson_config() {
  mps::ExperimentConfig config{
      .kind = mps::ExperimentKind::kShallowWater,
      .planet = {6371220, 7.292e-5, 9.80616, 287, 1004.5, 100000},
      .run = {0, 600, 60, 0},
      .grid = {6},
      .shallow_water = {.test_case = mps::ShallowWaterTestCase::kWilliamson5,
                        .scheme = mps::ShallowWaterScheme::kRusanov,
                        .reconstruction = mps::ReconstructionKind::kLinear,
                        .limiter = mps::LimiterKind::kBarthJespersen,
                        .cfl = .45,
                        .mean_depth_m = 5960,
                        .depth_floor_m = 10,
                        .flow_axis_z = 1,
                        .maximum_velocity_m_s = 20},
      .orography = {.kind = mps::OrographyKind::kWilliamson5},
      .diagnostics = {5},
      .output_directory = "x"};
  return config;
}

double terrain_force_for_bell(const double peak_height_m) {
  const mps::CubedSphereGrid grid(16, 6371220);
  constexpr double gravity = 9.80616;
  constexpr double rd = 287;
  constexpr double temperature = 280;
  mps::DryHydrostaticDerived derived{.cells = grid.cell_count(), .levels = 1};
  derived.pressure_pa.resize(derived.cells);
  derived.air_mass_kg_m2.assign(derived.cells, 1000);
  derived.velocity_m_s.assign(derived.cells, {});
  derived.temperature_k.assign(derived.cells, temperature);
  derived.geopotential_m2_s2.resize(derived.cells);
  for (std::size_t cell = 0; cell < derived.cells; ++cell) {
    const double phi = gravity * mps::linear_bell_surface_height_m(
                                     grid.cells()[cell].center, peak_height_m);
    derived.geopotential_m2_s2[cell] = phi;
    derived.pressure_pa[cell] = 100000 * std::exp(-phi / (rd * temperature));
  }
  const mps::PlanetParameters planet{grid.radius_m(), 0, gravity, rd, 1004.5, 100000};
  const auto source = mps::dry_hydrostatic_sources(grid, derived, planet);
  double maximum = 0;
  for (const auto value : source.geopotential_gradient_kg_m_s2)
    maximum = std::max(maximum, mps::norm(value));
  return maximum;
}
}  // namespace

MPS_TEST_CASE("Williamson 5 has a bounded short response and conserves mass") {
  const auto result = mps::run_shallow_water(williamson_config());
  MPS_CHECK(result.reached_end_time);
  MPS_CHECK(std::isfinite(result.final_diagnostics.energy));
  MPS_CHECK_NEAR(result.final_diagnostics.mass, result.initial_diagnostics.mass,
                 5e-13 * result.initial_diagnostics.mass);
  MPS_CHECK(result.final_diagnostics.energy != result.initial_diagnostics.energy);
}

MPS_TEST_CASE("hydrostatic bell terrain force is linear at small amplitude") {
  const double half = terrain_force_for_bell(5);
  const double full = terrain_force_for_bell(10);
  MPS_CHECK(half > 0);
  MPS_CHECK_NEAR(full / half, 2.0, 0.02);
}

int main() { return mps::test::run_all(); }
