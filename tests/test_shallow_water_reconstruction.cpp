#include <algorithm>
#include <cmath>
#include <vector>

#include "myplanetsim/dynamics/shallow_water_reconstruction.hpp"
#include "myplanetsim/dynamics/shallow_water_rhs.hpp"
#include "support/test.hpp"

namespace {

[[nodiscard]] mps::ShallowWaterState smooth_state(const mps::CubedSphereGrid& grid) {
  mps::ShallowWaterState state{.time_s = 0.0,
                               .step = 0,
                               .depth = std::vector<mps::Real>(grid.cell_count()),
                               .momentum = std::vector<mps::Vec3>(grid.cell_count())};
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
    const auto position = grid.cells()[cell].center;
    state.depth[cell] = 2.0 + 0.2 * position.x;
    const auto velocity =
        grid.radius_m() * mps::cross(mps::Vec3{0.1, -0.05, 0.2}, position);
    state.momentum[cell] = state.depth[cell] * velocity;
  }
  return state;
}

}  // namespace

MPS_TEST_CASE("linear shallow-water reconstruction preserves constant state") {
  const mps::CubedSphereGrid grid(8, 1.0);
  mps::ShallowWaterState state{.time_s = 0.0,
                               .step = 0,
                               .depth = std::vector<mps::Real>(grid.cell_count(), 2.0),
                               .momentum = std::vector<mps::Vec3>(grid.cell_count())};
  const auto faces = mps::reconstruct_shallow_water_face_states(
      grid, state, mps::ReconstructionKind::kLinear, mps::LimiterKind::kBarthJespersen,
      0.1);
  for (const auto& edge : faces.edges) {
    MPS_CHECK_EQ(edge.left.depth_m, 2.0);
    MPS_CHECK_EQ(edge.right.depth_m, 2.0);
    MPS_CHECK_EQ(mps::norm(edge.left.velocity_m_s), 0.0);
    MPS_CHECK_EQ(mps::norm(edge.right.velocity_m_s), 0.0);
  }
}

MPS_TEST_CASE("linear shallow-water limiter keeps near-dry faces above floor") {
  const mps::CubedSphereGrid grid(4, 1.0);
  auto state = smooth_state(grid);
  constexpr mps::Real floor = 0.5;
  state.depth[0] = std::nextafter(floor, 1.0);
  state.momentum[0] = {};
  const auto faces = mps::reconstruct_shallow_water_face_states(
      grid, state, mps::ReconstructionKind::kLinear, mps::LimiterKind::kBarthJespersen,
      floor);
  MPS_CHECK(faces.limiter_activations > 0);
  for (const auto& edge : faces.edges) {
    MPS_CHECK(edge.left.depth_m > floor);
    MPS_CHECK(edge.right.depth_m > floor);
  }
}

MPS_TEST_CASE("linear reconstructed RHS retains unique-edge mass cancellation") {
  const mps::CubedSphereGrid grid(8, 2.0);
  const auto state = smooth_state(grid);
  mps::ShallowWaterParameters parameters{};
  parameters.reconstruction = mps::ReconstructionKind::kLinear;
  parameters.limiter = mps::LimiterKind::kBarthJespersen;
  parameters.depth_floor_m = 0.1;
  const auto rhs =
      mps::assemble_shallow_water_rhs(grid, state, parameters, 3.0, {0.0, 0.0, 0.1});
  mps::Real sum = 0.0;
  mps::Real scale = 0.0;
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
    const mps::Real integrated = rhs.total.depth[cell] * grid.cells()[cell].area_m2;
    sum += integrated;
    scale += std::abs(integrated);
  }
  MPS_CHECK(std::abs(sum) < 5.0e-13 * scale);
}

int main() { return mps::test::run_all(); }
