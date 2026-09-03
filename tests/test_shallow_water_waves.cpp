#include <algorithm>
#include <cmath>
#include <vector>

#include "myplanetsim/diagnostics/reductions.hpp"
#include "myplanetsim/dynamics/shallow_water_benchmarks.hpp"
#include "myplanetsim/dynamics/shallow_water_driver.hpp"
#include "myplanetsim/numerics/explicit_steppers.hpp"
#include "support/test.hpp"

namespace {

[[nodiscard]] mps::ExperimentConfig wave_config(const mps::Index resolution) {
  mps::ExperimentConfig config{};
  config.kind = mps::ExperimentKind::kShallowWater;
  config.planet = {.radius_m = 1.0,
                   .rotation_rate_rad_s = 0.0,
                   .gravity_m_s2 = 1.0,
                   .gas_constant_j_kg_k = 1.0,
                   .heat_capacity_cp_j_kg_k = 2.0,
                   .reference_pressure_pa = 1.0};
  config.run = {
      .start_time_s = 0.0, .end_time_s = 0.1, .time_step_s = 0.002, .random_seed = 1};
  config.grid = {.cells_per_panel = resolution, .halo_width = 2};
  config.shallow_water.test_case = mps::ShallowWaterTestCase::kLinearWave;
  config.shallow_water.scheme = mps::ShallowWaterScheme::kRusanov;
  config.shallow_water.reconstruction = mps::ReconstructionKind::kLinear;
  config.shallow_water.limiter = mps::LimiterKind::kNone;
  config.shallow_water.cfl = 0.4;
  config.shallow_water.mean_depth_m = 1.0;
  config.shallow_water.depth_floor_m = 0.1;
  config.shallow_water.maximum_velocity_m_s = 0.0;
  config.output_directory = "output";
  return config;
}

[[nodiscard]] mps::Real wave_depth_error(const mps::Index resolution) {
  const auto config = wave_config(resolution);
  const auto result = mps::run_shallow_water(config);
  const mps::CubedSphereGrid grid(resolution, config.planet.radius_m);
  const auto exact = mps::make_linear_wave_state(grid, config.run.end_time_s,
                                                 config.planet.gravity_m_s2,
                                                 config.shallow_water.mean_depth_m);
  std::vector<mps::Real> weights(grid.cell_count());
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
    weights[cell] = grid.cells()[cell].area_m2;
  }
  return mps::diagnostics::weighted_error_norms(result.state.depth, exact.depth,
                                                weights)
      .l2;
}

[[nodiscard]] mps::Real inertial_error(const mps::Real time_step) {
  constexpr mps::Real coriolis = 1.3;
  constexpr mps::Real end_time = 1.0;
  std::vector<mps::Real> velocity{1.0, 0.0};
  mps::SspRk3 stepper(velocity.size());
  auto rhs = [=](const mps::Real, const std::span<const mps::Real> state,
                 const std::span<mps::Real> tendency) {
    tendency[0] = coriolis * state[1];
    tendency[1] = -coriolis * state[0];
  };
  mps::Real time = 0.0;
  while (time < end_time) {
    const mps::Real step = std::min(time_step, end_time - time);
    stepper.step(time, step, velocity, rhs);
    time += step;
  }
  return std::hypot(velocity[0] - std::cos(coriolis * end_time),
                    velocity[1] + std::sin(coriolis * end_time));
}

}  // namespace

MPS_TEST_CASE("nonrotating inertia-gravity wave error decreases with resolution") {
  const mps::Real coarse = wave_depth_error(8);
  const mps::Real fine = wave_depth_error(16);
  MPS_CHECK(fine < 0.8 * coarse);
}

MPS_TEST_CASE("Coriolis-only inertial oscillation has SSPRK3 time order") {
  const mps::Real coarse = inertial_error(0.05);
  const mps::Real fine = inertial_error(0.025);
  MPS_CHECK(fine < 0.14 * coarse);
}

MPS_TEST_CASE("geostrophic adjustment remains finite and depth-positive") {
  auto config = wave_config(8);
  config.planet.rotation_rate_rad_s = 0.2;
  config.shallow_water.test_case = mps::ShallowWaterTestCase::kGeostrophicAdjustment;
  const auto result = mps::run_shallow_water(config);
  MPS_CHECK(result.final_diagnostics.minimum_depth >
            config.shallow_water.depth_floor_m);
  MPS_CHECK(std::isfinite(result.final_diagnostics.energy));
}

MPS_TEST_CASE("uniform rest survives one thousand steps without new extrema") {
  auto config = wave_config(1);
  config.shallow_water.test_case = mps::ShallowWaterTestCase::kRest;
  config.run.end_time_s = 1.0;
  config.run.time_step_s = 0.001;
  const auto result = mps::run_shallow_water(config);
  MPS_CHECK_EQ(result.state.step, 1000U);
  MPS_CHECK_EQ(result.final_diagnostics.minimum_depth, 1.0);
  MPS_CHECK_EQ(result.final_diagnostics.maximum_depth, 1.0);
  MPS_CHECK_NEAR(result.final_diagnostics.mass, result.initial_diagnostics.mass,
                 5.0e-14 * result.initial_diagnostics.mass);
}

int main() { return mps::test::run_all(); }
