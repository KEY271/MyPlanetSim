#include <cmath>

#include "myplanetsim/dynamics/shallow_water_diffusion.hpp"
#include "myplanetsim/numerics/diffusion_stability.hpp"
#include "support/test.hpp"

namespace {

[[nodiscard]] mps::ShallowWaterState smooth_state(const mps::CubedSphereGrid& grid) {
  mps::ShallowWaterState state{.depth = std::vector<mps::Real>(grid.cell_count()),
                               .momentum = std::vector<mps::Vec3>(grid.cell_count())};
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
    const mps::Vec3 position = grid.cells()[cell].center;
    state.depth[cell] = 2.0 + 0.1 * position.x;
    const mps::Vec3 velocity = mps::cross(mps::Vec3{0.0, 0.0, 1.0}, position);
    state.momentum[cell] = state.depth[cell] * velocity;
  }
  return state;
}

}  // namespace

MPS_TEST_CASE("shallow-water diffusion is globally mass neutral") {
  const mps::CubedSphereGrid grid(12, 2.0);
  const auto state = smooth_state(grid);
  for (const auto kind :
       {mps::DiffusionKind::kLaplacian, mps::DiffusionKind::kBiharmonic}) {
    const auto tendency = mps::shallow_water_diffusion_tendency(grid, state, kind, 0.1);
    mps::Real mass_rate = 0.0;
    mps::Real absolute_rate = 0.0;
    for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
      const mps::Real integrated = grid.cells()[cell].area_m2 * tendency.depth[cell];
      mass_rate += integrated;
      absolute_rate += std::abs(integrated);
    }
    MPS_CHECK(std::abs(mass_rate) <= 5.0e-13 * absolute_rate);
  }
}

MPS_TEST_CASE("Laplacian and biharmonic diffusion damp a smooth depth mode") {
  const mps::CubedSphereGrid grid(12, 2.0);
  const auto state = smooth_state(grid);
  for (const auto kind :
       {mps::DiffusionKind::kLaplacian, mps::DiffusionKind::kBiharmonic}) {
    const auto tendency = mps::shallow_water_diffusion_tendency(grid, state, kind, 0.1);
    mps::Real variance_rate = 0.0;
    for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
      variance_rate +=
          grid.cells()[cell].area_m2 * (state.depth[cell] - 2.0) * tendency.depth[cell];
    }
    MPS_CHECK(variance_rate < 0.0);
  }
}

MPS_TEST_CASE("disabled diffusion is exactly zero") {
  const mps::CubedSphereGrid grid(4, 1.0);
  const auto tendency = mps::shallow_water_diffusion_tendency(
      grid, smooth_state(grid), mps::DiffusionKind::kNone, 0.0);
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
    MPS_CHECK_EQ(tendency.depth[cell], 0.0);
    MPS_CHECK_EQ(mps::norm_squared(tendency.momentum[cell]), 0.0);
  }
}

MPS_TEST_CASE("diffusion stability limit tightens under refinement") {
  const mps::CubedSphereGrid coarse(8, 1.0);
  const mps::CubedSphereGrid fine(16, 1.0);
  const auto stable_step = [](const mps::CubedSphereGrid& grid,
                              const mps::DiffusionKind kind) {
    return mps::stable_diffusion_time_step(grid, kind, 0.1, 10.0);
  };
  MPS_CHECK(stable_step(fine, mps::DiffusionKind::kLaplacian) <
            stable_step(coarse, mps::DiffusionKind::kLaplacian));
  MPS_CHECK(stable_step(fine, mps::DiffusionKind::kBiharmonic) <
            stable_step(coarse, mps::DiffusionKind::kBiharmonic));
  MPS_CHECK_EQ(
      mps::stable_diffusion_time_step(coarse, mps::DiffusionKind::kNone, 0.0, 10.0),
      10.0);
}

int main() { return mps::test::run_all(); }
