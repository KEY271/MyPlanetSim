#include <algorithm>
#include <cmath>
#include <vector>

#include "myplanetsim/dynamics/dry_hydrostatic_sources.hpp"
#include "myplanetsim/dynamics/surface_orography.hpp"
#include "support/test.hpp"

namespace {
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

MPS_TEST_CASE("hydrostatic bell terrain force is linear at small amplitude") {
  const double half = terrain_force_for_bell(5);
  const double full = terrain_force_for_bell(10);
  MPS_CHECK(half > 0);
  MPS_CHECK_NEAR(full / half, 2.0, 0.02);
}

int main() { return mps::test::run_all(); }
