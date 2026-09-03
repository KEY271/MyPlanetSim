#include <algorithm>
#include <bit>
#include <cmath>
#include <cstdint>
#include <vector>

#include "myplanetsim/diagnostics/reductions.hpp"
#include "myplanetsim/dynamics/shallow_water_benchmarks.hpp"
#include "myplanetsim/dynamics/shallow_water_compatible.hpp"
#include "myplanetsim/dynamics/shallow_water_driver.hpp"
#include "support/test.hpp"

namespace {

constexpr mps::Real kRadius = 6.37122e6;
constexpr mps::Real kGravity = 9.80616;
constexpr mps::Real kRotation = 7.292e-5;
constexpr mps::Real kDepth = 2.94e4 / kGravity;
constexpr mps::Real kVelocity =
    2.0 * 3.14159265358979323846 * kRadius / (12.0 * 86400.0);

[[nodiscard]] mps::ExperimentConfig config(const mps::Index resolution,
                                           const mps::ShallowWaterScheme scheme,
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
  value.shallow_water.test_case = mps::ShallowWaterTestCase::kWilliamson2;
  value.shallow_water.scheme = scheme;
  value.shallow_water.reconstruction = mps::ReconstructionKind::kLinear;
  value.shallow_water.limiter = mps::LimiterKind::kBarthJespersen;
  value.shallow_water.cfl = 0.45;
  value.shallow_water.mean_depth_m = kDepth;
  value.shallow_water.depth_floor_m = 1.0;
  value.shallow_water.maximum_velocity_m_s = kVelocity;
  value.output_directory = "output";
  return value;
}

[[nodiscard]] mps::Real momentum_median(const mps::CubedSphereGrid& grid,
                                        const mps::ShallowWaterTendency& tendency) {
  std::vector<mps::Real> magnitudes(grid.cell_count());
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
    magnitudes[cell] = mps::norm(tendency.momentum[cell]);
  }
  std::ranges::sort(magnitudes);
  return magnitudes[magnitudes.size() / 2];
}

[[nodiscard]] mps::Real drift(const mps::Real final_value,
                              const mps::Real initial_value) {
  return std::abs(
      mps::diagnostics::relative_drift(final_value, initial_value, initial_value));
}

[[nodiscard]] mps::Real williamson2_residual(const mps::Index resolution) {
  const mps::CubedSphereGrid grid(resolution, kRadius);
  const mps::CubedSphereDualTopology dual(grid);
  const auto state = mps::make_williamson2_state(grid, 0.0, kGravity, kRotation, kDepth,
                                                 kVelocity, {0.0, 0.0, 1.0});
  const auto rhs = mps::assemble_compatible_shallow_water_rhs(
      grid, dual, state,
      config(resolution, mps::ShallowWaterScheme::kCompatible, 600.0).shallow_water,
      kGravity, {0.0, 0.0, kRotation});
  return momentum_median(grid, rhs.total) / momentum_median(grid, rhs.pressure);
}

}  // namespace

MPS_TEST_CASE("compatible scheme keeps a rotating lake at rest") {
  const mps::CubedSphereGrid grid(8, kRadius);
  const mps::CubedSphereDualTopology dual(grid);
  const auto state = mps::make_resting_shallow_water_state(grid, 0.0, kDepth);
  mps::ShallowWaterParameters parameters{};
  parameters.scheme = mps::ShallowWaterScheme::kCompatible;
  parameters.depth_floor_m = 1.0;
  const auto rhs = mps::assemble_compatible_shallow_water_rhs(
      grid, dual, state, parameters, kGravity, {0.0, 0.0, kRotation});
  // Scale of a momentum tendency built from a free-surface gradient.
  const mps::Real momentum_scale = kGravity * kDepth * kDepth / kRadius;
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
    MPS_CHECK_NEAR(rhs.total.depth[cell], 0.0, 1.0e-18 * kDepth);
    MPS_CHECK_NEAR(mps::norm(rhs.total.momentum[cell]), 0.0, 1.0e-15 * momentum_scale);
  }
}

MPS_TEST_CASE("compatible mass flux cancels over the sphere") {
  const mps::CubedSphereGrid grid(8, kRadius);
  const mps::CubedSphereDualTopology dual(grid);
  const auto state =
      mps::make_williamson2_state(grid, 0.0, kGravity, kRotation, kDepth, kVelocity,
                                  mps::normalize(mps::Vec3{1.0, 1.0, 1.0}));
  mps::ShallowWaterParameters parameters{};
  parameters.scheme = mps::ShallowWaterScheme::kCompatible;
  parameters.depth_floor_m = 1.0;
  const auto rhs = mps::assemble_compatible_shallow_water_rhs(
      grid, dual, state, parameters, kGravity, {0.0, 0.0, kRotation});
  mps::Real tendency = 0.0;
  mps::Real scale = 0.0;
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
    const mps::Real weighted = grid.cells()[cell].area_m2 * rhs.flux.depth[cell];
    tendency += weighted;
    scale += std::abs(weighted);
  }
  MPS_CHECK(std::abs(tendency) < 5.0e-13 * scale);
}

MPS_TEST_CASE("compatible balance residual converges on Williamson 2") {
  const mps::Real coarse = williamson2_residual(12);
  const mps::Real fine = williamson2_residual(24);
  MPS_CHECK(coarse < 0.05);
  MPS_CHECK(fine < 0.5 * coarse);
}

MPS_TEST_CASE("compatible driver conserves mass and dissipates less than Rusanov") {
  const auto compatible =
      mps::run_shallow_water(config(12, mps::ShallowWaterScheme::kCompatible, 3600.0));
  const auto rusanov =
      mps::run_shallow_water(config(12, mps::ShallowWaterScheme::kRusanov, 3600.0));
  const mps::Real mass_drift = mps::diagnostics::relative_drift(
      compatible.final_diagnostics.mass, compatible.initial_diagnostics.mass,
      compatible.initial_diagnostics.mass);
  MPS_CHECK(std::abs(mass_drift) < 5.0e-13);
  MPS_CHECK(compatible.final_diagnostics.minimum_depth > 1.0);
  MPS_CHECK(
      drift(compatible.final_diagnostics.energy,
            compatible.initial_diagnostics.energy) <
      drift(rusanov.final_diagnostics.energy, rusanov.initial_diagnostics.energy));
  MPS_CHECK(drift(compatible.final_diagnostics.potential_enstrophy,
                  compatible.initial_diagnostics.potential_enstrophy) <
            drift(rusanov.final_diagnostics.potential_enstrophy,
                  rusanov.initial_diagnostics.potential_enstrophy));
}

MPS_TEST_CASE("compatible restart is bitwise continuous") {
  const auto value = config(8, mps::ShallowWaterScheme::kCompatible, 1800.0);
  const auto continuous = mps::run_shallow_water(value);
  const auto stopped = mps::run_shallow_water(value, std::nullopt, 2);
  MPS_CHECK(!stopped.reached_end_time);
  const auto restarted = mps::run_shallow_water(value, stopped.state);
  MPS_CHECK_EQ(restarted.state.step, continuous.state.step);
  for (std::size_t cell = 0; cell < continuous.state.cell_count(); ++cell) {
    MPS_CHECK_EQ(std::bit_cast<std::uint64_t>(restarted.state.depth[cell]),
                 std::bit_cast<std::uint64_t>(continuous.state.depth[cell]));
    for (std::size_t component = 0; component < 3; ++component) {
      MPS_CHECK_EQ(
          std::bit_cast<std::uint64_t>(restarted.state.momentum[cell][component]),
          std::bit_cast<std::uint64_t>(continuous.state.momentum[cell][component]));
    }
  }
}

int main() { return mps::test::run_all(); }
