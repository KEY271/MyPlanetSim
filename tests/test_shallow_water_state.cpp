#include <array>
#include <bit>
#include <cstdint>
#include <limits>
#include <stdexcept>

#include "myplanetsim/dynamics/shallow_water_state.hpp"
#include "support/test.hpp"

namespace {

[[nodiscard]] mps::ShallowWaterState sample_state(const mps::CubedSphereGrid& grid) {
  mps::ShallowWaterState state{.time_s = 2.0,
                               .step = 3,
                               .depth = std::vector<mps::Real>(grid.cell_count()),
                               .momentum = std::vector<mps::Vec3>(grid.cell_count())};
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
    state.depth[cell] = 10.0 + static_cast<mps::Real>(cell);
    state.momentum[cell] = mps::project_tangent(
        {1.0 + static_cast<mps::Real>(cell), -2.0, 3.0}, grid.cells()[cell].center);
  }
  return state;
}

}  // namespace

MPS_TEST_CASE("shallow-water state has a bitwise flat round trip") {
  const mps::CubedSphereGrid grid(2, 1.0);
  const auto original = sample_state(grid);
  const auto flat = mps::flatten_shallow_water_state(original);
  const auto restored = mps::unflatten_shallow_water_state(
      original.time_s, original.step, flat, grid.cell_count());
  MPS_CHECK_EQ(restored.step, original.step);
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
    MPS_CHECK_EQ(std::bit_cast<std::uint64_t>(restored.depth[cell]),
                 std::bit_cast<std::uint64_t>(original.depth[cell]));
    for (std::size_t component = 0; component < 3; ++component) {
      MPS_CHECK_EQ(std::bit_cast<std::uint64_t>(restored.momentum[cell][component]),
                   std::bit_cast<std::uint64_t>(original.momentum[cell][component]));
    }
  }
  mps::validate_shallow_water_state(grid, restored, 0.1);
}

MPS_TEST_CASE("shallow-water state rejects invalid shape depth and momentum") {
  const mps::CubedSphereGrid grid(2, 1.0);
  auto state = sample_state(grid);
  state.depth.pop_back();
  MPS_CHECK_THROWS_AS(mps::validate_shallow_water_state(grid, state, 0.1),
                      std::invalid_argument);
  state = sample_state(grid);
  state.depth[0] = 0.1;
  MPS_CHECK_THROWS_AS(mps::validate_shallow_water_state(grid, state, 0.1),
                      std::runtime_error);
  state = sample_state(grid);
  state.momentum[0] = grid.cells()[0].center;
  MPS_CHECK_THROWS_AS(mps::validate_shallow_water_state(grid, state, 0.1),
                      std::runtime_error);
  state = sample_state(grid);
  state.momentum[0].x = std::numeric_limits<mps::Real>::quiet_NaN();
  MPS_CHECK_THROWS_AS(mps::validate_shallow_water_state(grid, state, 0.1),
                      std::runtime_error);
}

MPS_TEST_CASE("shallow-water momentum projection restores tangency") {
  const mps::CubedSphereGrid grid(1, 1.0);
  auto state = sample_state(grid);
  state.momentum[0] = state.momentum[0] + grid.cells()[0].center;
  mps::project_shallow_water_momentum(grid, state);
  mps::validate_shallow_water_state(grid, state, 0.1);
}

MPS_TEST_CASE("shallow-water flat state rejects layout shape mismatch") {
  const std::array<mps::Real, 2> invalid{1.0, 2.0};
  MPS_CHECK_THROWS_AS(mps::unflatten_shallow_water_state(0.0, 0, invalid, 1),
                      std::invalid_argument);
}

int main() { return mps::test::run_all(); }
