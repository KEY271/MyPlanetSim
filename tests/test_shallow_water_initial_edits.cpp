#include <cmath>
#include <numbers>
#include <stdexcept>
#include <vector>

#include "myplanetsim/dynamics/shallow_water_benchmarks.hpp"
#include "myplanetsim/dynamics/shallow_water_initial_edits.hpp"
#include "support/test.hpp"

namespace {

[[nodiscard]] mps::ExperimentConfig config() {
  mps::ExperimentConfig value{};
  value.kind = mps::ExperimentKind::kShallowWater;
  value.planet = {.radius_m = 1.0, .rotation_rate_rad_s = 0.0, .gravity_m_s2 = 1.0,
                  .gas_constant_j_kg_k = 1.0, .heat_capacity_cp_j_kg_k = 2.0,
                  .reference_pressure_pa = 1.0};
  value.run = {.start_time_s = 0.0, .end_time_s = 1.0, .time_step_s = 0.1,
               .random_seed = 1};
  value.grid = {.cells_per_panel = 4};
  value.shallow_water.test_case = mps::ShallowWaterTestCase::kRest;
  value.shallow_water.cfl = 0.5;
  value.shallow_water.mean_depth_m = 10.0;
  value.shallow_water.depth_floor_m = 0.1;
  value.output_directory = "output";
  return value;
}

}  // namespace

MPS_TEST_CASE("preserved Gaussian edit has zero global mass change") {
  const auto settings = config();
  const mps::CubedSphereGrid grid(4, settings.planet.radius_m);
  auto state = mps::make_shallow_water_initial_state(grid, settings);
  const auto edit = mps::InitialConditionEditV1{
      .center_unit = mps::normalize({1.0, 2.0, 3.0}),
      .amplitude_m = 1.0,
      .sigma_rad = 0.2,
      .mass_policy = mps::MassPolicy::kPreserveGlobal};
  const auto diagnostics = mps::apply_initial_condition_edits(grid, state, {&edit, 1}, 0.1);
  MPS_CHECK_NEAR(diagnostics.added_volume_m3, 0.0, 2.0e-13 * grid.total_area_m2());
  MPS_CHECK(diagnostics.maximum_depth_m > settings.shallow_water.mean_depth_m);
  MPS_CHECK(diagnostics.minimum_depth_m < settings.shallow_water.mean_depth_m);
}

MPS_TEST_CASE("allow-change Gaussian edit records added volume") {
  const auto settings = config();
  const mps::CubedSphereGrid grid(4, settings.planet.radius_m);
  auto state = mps::make_shallow_water_initial_state(grid, settings);
  const auto before = state.depth;
  const auto edit = mps::InitialConditionEditV1{
      .center_unit = {0.0, 0.0, 1.0}, .amplitude_m = 0.5, .sigma_rad = 0.3,
      .mass_policy = mps::MassPolicy::kAllowChange};
  const auto diagnostics = mps::apply_initial_condition_edits(grid, state, {&edit, 1}, 0.1);
  MPS_CHECK(diagnostics.added_volume_m3 > 0.0);
  MPS_CHECK(state.depth != before);
}

MPS_TEST_CASE("Gaussian edit keeps momentum tangent and rejects invalid edit") {
  const auto settings = config();
  const mps::CubedSphereGrid grid(2, settings.planet.radius_m);
  auto state = mps::make_shallow_water_initial_state(grid, settings);
  const auto invalid = mps::InitialConditionEditV1{
      .center_unit = {0.0, 0.0, 2.0}, .amplitude_m = 1.0, .sigma_rad = 0.2};
  MPS_CHECK_THROWS_AS(mps::apply_initial_condition_edits(grid, state, {&invalid, 1}, 0.1),
                      std::invalid_argument);
  const auto edit = mps::InitialConditionEditV1{
      .center_unit = grid.cells()[0].center, .amplitude_m = 0.01, .sigma_rad = 0.4};
  const auto diagnostics = mps::apply_initial_condition_edits(grid, state, {&edit, 1}, 0.1);
  (void)diagnostics;
  mps::validate_shallow_water_state(grid, state, 0.1);
}

MPS_TEST_CASE("multiple edits are applied in order") {
  const auto settings = config();
  const mps::CubedSphereGrid grid(2, settings.planet.radius_m);
  auto state = mps::make_shallow_water_initial_state(grid, settings);
  const std::vector<mps::InitialConditionEditV1> edits{
      {.center_unit = grid.cells()[0].center, .amplitude_m = 0.1, .sigma_rad = 0.3},
      {.center_unit = grid.cells()[1].center, .amplitude_m = 0.2, .sigma_rad = 0.3}};
  const auto diagnostics = mps::apply_initial_condition_edits(grid, state, edits, 0.1);
  (void)diagnostics;
  MPS_CHECK(std::isfinite(state.depth[0]));
}

int main() { return mps::test::run_all(); }
