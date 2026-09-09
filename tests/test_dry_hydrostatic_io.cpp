#include <sstream>
#include <stdexcept>

#include "myplanetsim/diagnostics/dry_hydrostatic_diagnostics.hpp"
#include "myplanetsim/io/checkpoint.hpp"
#include "myplanetsim/io/run_metadata.hpp"
#include "support/test.hpp"

namespace {

[[nodiscard]] mps::ExperimentConfig uniform_dry_config(const mps::Index levels) {
  mps::ExperimentConfig config;
  config.kind = mps::ExperimentKind::kDryHydrostatic;
  config.planet.radius_m = 6371220.0;
  config.planet.rotation_rate_rad_s = 7.29212e-5;
  config.planet.gravity_m_s2 = 9.80616;
  config.planet.gas_constant_j_kg_k = 287.0;
  config.planet.heat_capacity_cp_j_kg_k = 1004.0;
  config.planet.reference_pressure_pa = 100000.0;
  config.grid.cells_per_panel = 2;
  config.run.end_time_s = 60.0;
  config.run.time_step_s = 10.0;
  const auto coefficients = mps::uniform_sigma_coefficients(1000.0, levels);
  config.vertical.levels = levels;
  config.vertical.a_half_pa = coefficients.a_half_pa;
  config.vertical.b_half = coefficients.b_half;
  config.vertical.surface_pressure_pa = 100000.0;
  config.vertical.minimum_surface_pressure_pa = 80000.0;
  config.vertical.maximum_surface_pressure_pa = 120000.0;
  config.vertical.minimum_pressure_thickness_pa = 100.0;
  config.vertical.initial_temperature_k = 288.0;
  config.vertical.initial_potential_temperature_k = 300.0;
  config.vertical.temperature_floor_k = 100.0;
  config.output_directory = "output";
  return config;
}

}  // namespace
MPS_TEST_CASE("dry state checkpoint layout and budgets are deterministic") {
  mps::DryHydrostaticState s{.surface_pressure_pa = {100000},
                             .horizontal_momentum_mass_kg_m_s = {{0, 0, 0}},
                             .potential_temperature_mass_k_kg_m2 = {3000},
                             .tracer_mass_kg_m2 = {2}};
  auto flat = mps::flatten_dry_hydrostatic_state(s, 1);
  std::stringstream io;
  mps::write_checkpoint(io, {0, 0, flat, "fingerprint",
                             std::string(mps::kDryHydrostaticCheckpointLayout)});
  auto c = mps::read_checkpoint(io, "fingerprint", mps::kDryHydrostaticCheckpointLayout,
                                flat.size());
  MPS_CHECK_EQ(c.state.size(), flat.size());
  mps::CubedSphereGrid g(1, 2);
}
MPS_TEST_CASE("terrain changes checkpoint identity without changing its layout") {
  const auto flat = uniform_dry_config(2);
  auto terrain = flat;
  terrain.dry_hydrostatic.test_case = mps::DryHydrostaticTestCase::kLinearMountainWave;
  terrain.orography.kind = mps::OrographyKind::kLinearBell;
  const auto flat_fingerprint = mps::config_fingerprint(flat);
  const auto terrain_fingerprint = mps::config_fingerprint(terrain);
  MPS_CHECK(flat_fingerprint != terrain_fingerprint);
  std::stringstream stream;
  mps::write_checkpoint(
      stream, {.time_s = 0,
               .step = 0,
               .state = {1},
               .config_fingerprint = flat_fingerprint,
               .layout_id = std::string(mps::kDryHydrostaticCheckpointLayout)});
  MPS_CHECK_THROWS_AS(mps::read_checkpoint(stream, terrain_fingerprint,
                                           mps::kDryHydrostaticCheckpointLayout, 1),
                      std::runtime_error);
}
int main() { return mps::test::run_all(); }
