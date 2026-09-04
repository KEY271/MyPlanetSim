#include <bit>
#include <cstdint>
#include <sstream>
#include <vector>

#include "myplanetsim/dynamics/shallow_water_driver.hpp"
#include "support/test.hpp"

namespace {

[[nodiscard]] mps::ExperimentConfig config() {
  mps::ExperimentConfig value{};
  value.kind = mps::ExperimentKind::kShallowWater;
  value.planet = {.radius_m = 2.0,
                  .rotation_rate_rad_s = 0.1,
                  .gravity_m_s2 = 3.0,
                  .gas_constant_j_kg_k = 4.0,
                  .heat_capacity_cp_j_kg_k = 5.0,
                  .reference_pressure_pa = 6.0};
  value.run = {
      .start_time_s = 0.0, .end_time_s = 1.0, .time_step_s = 0.2, .random_seed = 7};
  value.grid = {.cells_per_panel = 4};
  value.shallow_water.test_case = mps::ShallowWaterTestCase::kRest;
  value.shallow_water.scheme = mps::ShallowWaterScheme::kRusanov;
  value.shallow_water.reconstruction = mps::ReconstructionKind::kPiecewiseConstant;
  value.shallow_water.limiter = mps::LimiterKind::kNone;
  value.shallow_water.cfl = 0.5;
  value.shallow_water.mean_depth_m = 2.0;
  value.shallow_water.depth_floor_m = 0.1;
  value.shallow_water.maximum_velocity_m_s = 0.0;
  value.output_directory = "output";
  return value;
}

}  // namespace

struct FrameTrace {
  std::vector<std::uint64_t> steps;
  std::uint64_t cancel_after = 0;
};

void record_frame(const mps::ShallowWaterState& state, void* context) {
  auto& trace = *static_cast<FrameTrace*>(context);
  trace.steps.push_back(state.step);
}

bool cancel_at_step(void* context) {
  const auto& trace = *static_cast<const FrameTrace*>(context);
  return trace.steps.size() > 1 && trace.steps.back() >= trace.cancel_after;
}

void throw_from_frame(const mps::ShallowWaterState&, void*) {
  throw std::runtime_error("observer failure");
}

MPS_TEST_CASE("gravity-wave CFL respects configured bounds") {
  const mps::CubedSphereGrid grid(4, 2.0);
  const mps::ShallowWaterState state{
      .depth = std::vector<mps::Real>(grid.cell_count(), 2.0),
      .momentum = std::vector<mps::Vec3>(grid.cell_count())};
  const mps::Real half =
      mps::stable_shallow_water_time_step(grid, state, 3.0, 0.5, 10.0);
  const mps::Real full =
      mps::stable_shallow_water_time_step(grid, state, 3.0, 1.0, 10.0);
  MPS_CHECK_NEAR(full, 2.0 * half, 2.0e-15);
  MPS_CHECK_NEAR(mps::shallow_water_cfl_number(grid, state, 3.0, half), 0.5, 2.0e-15);
  MPS_CHECK(mps::stable_shallow_water_time_step(grid, state, 3.0, 0.5, 0.01) <= 0.01);
}

MPS_TEST_CASE("shallow-water driver lands exactly on final time") {
  const auto result = mps::run_shallow_water(config());
  MPS_CHECK(result.reached_end_time);
  MPS_CHECK_EQ(result.state.time_s, config().run.end_time_s);
  MPS_CHECK(result.maximum_cfl <= config().shallow_water.cfl + 5.0e-15);
  MPS_CHECK_NEAR(result.final_diagnostics.mass, result.initial_diagnostics.mass,
                 5.0e-14 * result.initial_diagnostics.mass);
}

MPS_TEST_CASE("shallow-water stopped restart is bitwise continuous") {
  const auto continuous = mps::run_shallow_water(config());
  const auto stopped = mps::run_shallow_water(config(), std::nullopt, 2);
  MPS_CHECK(!stopped.reached_end_time);
  const auto restarted = mps::run_shallow_water(config(), stopped.state);
  MPS_CHECK_EQ(restarted.state.step, continuous.state.step);
  MPS_CHECK_EQ(std::bit_cast<std::uint64_t>(restarted.state.time_s),
               std::bit_cast<std::uint64_t>(continuous.state.time_s));
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

MPS_TEST_CASE("shallow-water diagnostics are sampled on the configured interval") {
  auto value = config();
  value.diagnostics.interval_steps = 2;
  const auto result = mps::run_shallow_water(value);
  MPS_CHECK(result.samples.size() >= 2);
  MPS_CHECK_EQ(result.samples.front().step, std::uint64_t{0});
  MPS_CHECK_EQ(result.samples.back().step, result.state.step);
  MPS_CHECK_EQ(result.samples.back().time_s, result.state.time_s);
  for (std::size_t index = 1; index + 1 < result.samples.size(); ++index) {
    MPS_CHECK_EQ(result.samples[index].step % 2, std::uint64_t{0});
  }
  MPS_CHECK_EQ(result.samples.back().invariants.mass, result.final_diagnostics.mass);
}

MPS_TEST_CASE("edited initial diagnostics describe the pre-integration state") {
  const std::vector<mps::InitialConditionEditV1> edits{{
      .center_unit = {1.0, 0.0, 0.0},
      .amplitude_m = 0.1,
      .sigma_rad = 0.3,
      .mass_policy = mps::MassPolicy::kPreserveGlobal,
  }};
  const auto result =
      mps::run_shallow_water(config(), std::nullopt, std::nullopt, edits);
  MPS_CHECK_EQ(result.initial_diagnostics.mass, result.samples.front().invariants.mass);
  MPS_CHECK_EQ(result.initial_diagnostics.energy,
               result.samples.front().invariants.energy);
  MPS_CHECK_EQ(result.initial_diagnostics.potential_enstrophy,
               result.samples.front().invariants.potential_enstrophy);
}

MPS_TEST_CASE("frame observer reports initial intervals and one final frame") {
  auto value = config();
  value.run.end_time_s = 1.1;
  value.diagnostics.interval_steps = 2;
  FrameTrace trace;
  const auto result =
      mps::run_shallow_water(value, std::nullopt, std::nullopt, {},
                             {.on_frame = record_frame, .observer_context = &trace});
  MPS_CHECK_EQ(trace.steps.front(), std::uint64_t{0});
  MPS_CHECK_EQ(trace.steps.back(), result.state.step);
  for (std::size_t index = 1; index < trace.steps.size(); ++index) {
    MPS_CHECK(trace.steps[index] > trace.steps[index - 1]);
  }
  for (std::size_t index = 1; index + 1 < trace.steps.size(); ++index) {
    MPS_CHECK_EQ(trace.steps[index] % 2, std::uint64_t{0});
  }
}

MPS_TEST_CASE("cancellation is evaluated at step boundaries") {
  auto value = config();
  value.run.end_time_s = 10.0;
  FrameTrace trace{.steps = {}, .cancel_after = 2};
  const auto result = mps::run_shallow_water(value, std::nullopt, std::nullopt, {},
                                             {.on_frame = record_frame,
                                              .observer_context = &trace,
                                              .is_cancelled = cancel_at_step,
                                              .cancellation_context = &trace});
  MPS_CHECK(!result.reached_end_time);
  MPS_CHECK_EQ(result.state.step, std::uint64_t{2});
  MPS_CHECK_EQ(trace.steps.back(), result.state.step);
}

MPS_TEST_CASE("observer exceptions are propagated") {
  MPS_CHECK_THROWS_AS(mps::run_shallow_water(config(), std::nullopt, std::nullopt, {},
                                             {.on_frame = throw_from_frame}),
                      std::runtime_error);
}

int main() { return mps::test::run_all(); }
