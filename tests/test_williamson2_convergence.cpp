#include <cmath>
#include <vector>

#include "myplanetsim/diagnostics/reductions.hpp"
#include "myplanetsim/dynamics/shallow_water_benchmarks.hpp"
#include "myplanetsim/dynamics/shallow_water_driver.hpp"
#include "support/test.hpp"

namespace {

struct WilliamsonErrors {
  mps::Real depth_l2;
  mps::Real velocity_l2;
  mps::Real relative_mass_drift;
};

[[nodiscard]] mps::ExperimentConfig williamson_config(const mps::Index resolution,
                                                      const mps::Vec3 axis) {
  mps::ExperimentConfig config{};
  config.kind = mps::ExperimentKind::kShallowWater;
  config.planet = {.radius_m = 6.37122e6,
                   .rotation_rate_rad_s = 7.292e-5,
                   .gravity_m_s2 = 9.80616,
                   .gas_constant_j_kg_k = 287.0,
                   .heat_capacity_cp_j_kg_k = 1004.5,
                   .reference_pressure_pa = 1.0e5};
  config.run = {.start_time_s = 0.0,
                .end_time_s = 3600.0,
                .time_step_s = 300.0,
                .random_seed = 1};
  config.grid = {.cells_per_panel = resolution, .halo_width = 2};
  config.shallow_water.test_case = mps::ShallowWaterTestCase::kWilliamson2;
  config.shallow_water.scheme = mps::ShallowWaterScheme::kRusanov;
  config.shallow_water.reconstruction = mps::ReconstructionKind::kLinear;
  config.shallow_water.limiter = mps::LimiterKind::kBarthJespersen;
  config.shallow_water.cfl = 0.45;
  config.shallow_water.mean_depth_m = 2.94e4 / config.planet.gravity_m_s2;
  config.shallow_water.depth_floor_m = 1.0;
  config.shallow_water.flow_axis_x = axis.x;
  config.shallow_water.flow_axis_y = axis.y;
  config.shallow_water.flow_axis_z = axis.z;
  config.shallow_water.maximum_velocity_m_s =
      2.0 * 3.14159265358979323846 * config.planet.radius_m / (12.0 * 86400.0);
  config.output_directory = "output";
  return config;
}

[[nodiscard]] WilliamsonErrors errors(const mps::Index resolution,
                                      const mps::Vec3 axis) {
  const auto config = williamson_config(resolution, axis);
  const auto result = mps::run_shallow_water(config);
  const mps::CubedSphereGrid grid(resolution, config.planet.radius_m);
  const auto exact = mps::make_williamson2_state(
      grid, config.run.end_time_s, config.planet.gravity_m_s2,
      config.planet.rotation_rate_rad_s, config.shallow_water.mean_depth_m,
      config.shallow_water.maximum_velocity_m_s, axis);
  std::vector<mps::Real> weights(grid.cell_count());
  std::vector<mps::Real> velocity_error(grid.cell_count());
  std::vector<mps::Real> zero(grid.cell_count());
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
    weights[cell] = grid.cells()[cell].area_m2;
    velocity_error[cell] =
        mps::norm(result.state.velocity(cell) - exact.velocity(cell));
  }
  return {
      .depth_l2 = mps::diagnostics::weighted_error_norms(result.state.depth,
                                                         exact.depth, weights)
                      .l2,
      .velocity_l2 =
          mps::diagnostics::weighted_error_norms(velocity_error, zero, weights).l2,
      .relative_mass_drift = mps::diagnostics::relative_drift(
          result.final_diagnostics.mass, result.initial_diagnostics.mass,
          result.initial_diagnostics.mass),
  };
}

}  // namespace

MPS_TEST_CASE("Williamson 2 short-run error decreases under refinement") {
  const auto coarse = errors(6, {0.0, 0.0, 1.0});
  const auto fine = errors(12, {0.0, 0.0, 1.0});
  MPS_CHECK(fine.depth_l2 < 0.8 * coarse.depth_l2);
  MPS_CHECK(fine.velocity_l2 < 0.8 * coarse.velocity_l2);
  MPS_CHECK(std::abs(coarse.relative_mass_drift) < 5.0e-13);
  MPS_CHECK(std::abs(fine.relative_mass_drift) < 5.0e-13);
}

MPS_TEST_CASE("Williamson 2 oblique-axis error remains comparable") {
  const auto aligned = errors(8, {0.0, 0.0, 1.0});
  const auto oblique = errors(8, mps::normalize(mps::Vec3{1.0, 1.0, 1.0}));
  MPS_CHECK(oblique.depth_l2 < 2.0 * aligned.depth_l2);
  MPS_CHECK(oblique.velocity_l2 < 2.0 * aligned.velocity_l2);
  MPS_CHECK(std::abs(oblique.relative_mass_drift) < 5.0e-13);
}

int main() { return mps::test::run_all(); }
