#include <algorithm>
#include <cmath>

#include "myplanetsim/dynamics/shallow_water_benchmarks.hpp"
#include "myplanetsim/dynamics/shallow_water_driver.hpp"
#include "myplanetsim/dynamics/shallow_water_rhs.hpp"
#include "support/test.hpp"

namespace {

constexpr mps::Real kRadius = 6.37122e6;
constexpr mps::Real kGravity = 9.80616;
constexpr mps::Real kRotation = 7.292e-5;
constexpr mps::Real kDepth = 10000.0;
constexpr mps::Real kVelocity = 80.0;
constexpr mps::Vec3 kAxis{0.0, 0.0, 1.0};

[[nodiscard]] mps::ExperimentConfig config(const mps::Index resolution,
                                           const mps::Real end_time_s) {
  mps::ExperimentConfig value{};
  value.kind = mps::ExperimentKind::kShallowWater;
  value.planet = {.radius_m = kRadius,
                  .rotation_rate_rad_s = kRotation,
                  .gravity_m_s2 = kGravity,
                  .gas_constant_j_kg_k = 287.0,
                  .heat_capacity_cp_j_kg_k = 1004.5,
                  .reference_pressure_pa = 1.0e5};
  value.run = {.start_time_s = 0.0,
               .end_time_s = end_time_s,
               .time_step_s = 300.0,
               .random_seed = 1};
  value.grid = {.cells_per_panel = resolution, .halo_width = 2};
  value.shallow_water.test_case = mps::ShallowWaterTestCase::kGalewsky;
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

[[nodiscard]] mps::ShallowWaterState balanced_state(const mps::CubedSphereGrid& grid) {
  return mps::make_galewsky_state(grid, 0.0, kGravity, kRotation, kDepth, kVelocity,
                                  kAxis, false);
}

// Area-weighted RMS of a momentum tendency, used as a balance-residual measure.
[[nodiscard]] mps::Real momentum_rms(const mps::CubedSphereGrid& grid,
                                     const mps::ShallowWaterTendency& tendency) {
  mps::Real weighted = 0.0;
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
    weighted += grid.cells()[cell].area_m2 * mps::norm_squared(tendency.momentum[cell]);
  }
  return std::sqrt(weighted / grid.total_area_m2());
}

[[nodiscard]] mps::Real balance_residual_ratio(const mps::Index resolution) {
  const mps::CubedSphereGrid grid(resolution, kRadius);
  const auto state = balanced_state(grid);
  const auto rhs = mps::assemble_shallow_water_rhs(
      grid, state, config(resolution, 600.0).shallow_water, kGravity, kRotation);
  return momentum_rms(grid, rhs.total) / momentum_rms(grid, rhs.coriolis);
}

}  // namespace

MPS_TEST_CASE(
    "Galewsky jet is confined, positive and matches the requested mean depth") {
  const mps::CubedSphereGrid grid(24, kRadius);
  const auto state = balanced_state(grid);
  mps::validate_shallow_water_state(grid, state, 1.0, kDepth * kVelocity);
  mps::Real mass = 0.0;
  mps::Real maximum_speed = 0.0;
  mps::Real southern_speed = 0.0;
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
    const mps::Real speed = mps::norm(state.velocity(cell));
    mass += grid.cells()[cell].area_m2 * state.depth[cell];
    maximum_speed = std::max(maximum_speed, speed);
    if (grid.cells()[cell].center.z <= 0.0) {
      southern_speed = std::max(southern_speed, speed);
    }
  }
  MPS_CHECK_NEAR(mass / grid.total_area_m2(), kDepth, 1.0e-9 * kDepth);
  MPS_CHECK(maximum_speed > 0.9 * kVelocity);
  MPS_CHECK(maximum_speed <= kVelocity);
  MPS_CHECK_EQ(southern_speed, 0.0);
}

MPS_TEST_CASE("Galewsky height perturbation is local and leaves the jet velocity") {
  const mps::CubedSphereGrid grid(24, kRadius);
  const auto balanced = balanced_state(grid);
  const auto perturbed = mps::make_galewsky_state(grid, 0.0, kGravity, kRotation,
                                                  kDepth, kVelocity, kAxis, true);
  mps::Real maximum_difference = 0.0;
  mps::Real southern_difference = 0.0;
  mps::Real maximum_velocity_difference = 0.0;
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
    const mps::Real difference = perturbed.depth[cell] - balanced.depth[cell];
    maximum_difference = std::max(maximum_difference, std::abs(difference));
    if (grid.cells()[cell].center.z <= 0.0) {
      southern_difference = std::max(southern_difference, std::abs(difference));
    }
    maximum_velocity_difference =
        std::max(maximum_velocity_difference,
                 mps::norm(perturbed.velocity(cell) - balanced.velocity(cell)));
  }
  MPS_CHECK(maximum_difference > 40.0);
  MPS_CHECK(maximum_difference < 120.0);
  MPS_CHECK(southern_difference < 1.0);
  MPS_CHECK(maximum_velocity_difference < 1.0e-9 * kVelocity);
}

MPS_TEST_CASE("Galewsky initial balance residual decreases under refinement") {
  const mps::Real coarse = balance_residual_ratio(24);
  const mps::Real fine = balance_residual_ratio(48);
  MPS_CHECK(coarse < 1.0);
  MPS_CHECK(fine < 0.4 * coarse);
}

MPS_TEST_CASE("Galewsky short integration conserves mass and keeps depth positive") {
  const auto value = config(8, 3600.0);
  const auto result = mps::run_shallow_water(value);
  MPS_CHECK(result.reached_end_time);
  MPS_CHECK_NEAR(result.final_diagnostics.mass, result.initial_diagnostics.mass,
                 5.0e-13 * result.initial_diagnostics.mass);
  MPS_CHECK(result.final_diagnostics.minimum_depth > value.shallow_water.depth_floor_m);
}

MPS_TEST_CASE("Galewsky balanced jet does not imprint the cubed-sphere wave four") {
  const auto value = config(8, 3600.0);
  const mps::CubedSphereGrid grid(value.grid.cells_per_panel, kRadius);
  const auto result = mps::run_shallow_water(value, balanced_state(grid));
  const auto initial =
      mps::diagnose_depth_wave_mode(grid, balanced_state(grid), kAxis, 4);
  const auto final_wave = mps::diagnose_depth_wave_mode(grid, result.state, kAxis, 4);
  MPS_CHECK(final_wave.amplitude < 1.5 * initial.amplitude);
}

int main() { return mps::test::run_all(); }
