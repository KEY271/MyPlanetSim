#include <cmath>

#include "myplanetsim/diagnostics/shallow_water_diagnostics.hpp"
#include "myplanetsim/dynamics/shallow_water_benchmarks.hpp"
#include "myplanetsim/dynamics/shallow_water_driver.hpp"
#include "support/test.hpp"

namespace {

constexpr mps::Real kRadius = 6.37122e6;
constexpr mps::Real kGravity = 9.80616;
constexpr mps::Real kRotation = 7.292e-5;
constexpr mps::Real kDepth = 8000.0;
constexpr mps::Real kVelocity = 50.0;

[[nodiscard]] mps::ExperimentConfig config(const mps::Index resolution) {
  mps::ExperimentConfig value{};
  value.kind = mps::ExperimentKind::kShallowWater;
  value.planet = {.radius_m = kRadius,
                  .rotation_rate_rad_s = kRotation,
                  .gravity_m_s2 = kGravity,
                  .gas_constant_j_kg_k = 287.0,
                  .heat_capacity_cp_j_kg_k = 1004.5,
                  .reference_pressure_pa = 1.0e5};
  value.run = {
      .start_time_s = 0.0, .end_time_s = 600.0, .time_step_s = 300.0, .random_seed = 1};
  value.grid = {.cells_per_panel = resolution};
  value.shallow_water.test_case = mps::ShallowWaterTestCase::kWilliamson6;
  value.shallow_water.scheme = mps::ShallowWaterScheme::kRusanov;
  value.shallow_water.reconstruction = mps::ReconstructionKind::kLinear;
  value.shallow_water.limiter = mps::LimiterKind::kBarthJespersen;
  value.shallow_water.cfl = 0.45;
  value.shallow_water.mean_depth_m = kDepth;
  value.shallow_water.depth_floor_m = 1.0;
  value.shallow_water.maximum_velocity_m_s = kVelocity;
  value.output_directory = "output";
  return value;
}

}  // namespace

MPS_TEST_CASE("Williamson 6 initial state is positive finite and tangent") {
  const mps::CubedSphereGrid grid(16, kRadius);
  const auto state = mps::make_williamson6_state(grid, 0.0, kGravity, kRotation, kDepth,
                                                 kVelocity, {0.0, 0.0, 1.0});
  mps::validate_shallow_water_state(grid, state, 1.0, kDepth * kVelocity);
  const auto wave4 = mps::diagnose_depth_wave_mode(grid, state, {0.0, 0.0, 1.0}, 4);
  const auto wave3 = mps::diagnose_depth_wave_mode(grid, state, {0.0, 0.0, 1.0}, 3);
  MPS_CHECK(wave4.amplitude > 100.0 * wave3.amplitude);
  MPS_CHECK(std::abs(wave4.phase_rad) < 1.0e-12);
}

MPS_TEST_CASE("Williamson 6 short integration preserves mass and wave four") {
  const auto value = config(8);
  const mps::CubedSphereGrid grid(value.grid.cells_per_panel, kRadius);
  const auto initial = mps::make_shallow_water_initial_state(grid, value);
  const auto initial_wave =
      mps::diagnose_depth_wave_mode(grid, initial, {0.0, 0.0, 1.0}, 4);
  const auto result = mps::run_shallow_water(value);
  const auto final_wave =
      mps::diagnose_depth_wave_mode(grid, result.state, {0.0, 0.0, 1.0}, 4);
  MPS_CHECK_NEAR(result.final_diagnostics.mass, result.initial_diagnostics.mass,
                 5.0e-13 * result.initial_diagnostics.mass);
  MPS_CHECK(final_wave.amplitude > 0.9 * initial_wave.amplitude);
  MPS_CHECK(result.final_diagnostics.minimum_depth > 1.0);
}

int main() { return mps::test::run_all(); }
