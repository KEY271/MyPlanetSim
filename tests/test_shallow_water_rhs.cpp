#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

#include "myplanetsim/diagnostics/reductions.hpp"
#include "myplanetsim/dynamics/shallow_water_rhs.hpp"
#include "support/test.hpp"

namespace {

[[nodiscard]] mps::ShallowWaterState make_state(const mps::CubedSphereGrid& grid,
                                                const mps::Real depth) {
  return {.time_s = 0.0,
          .step = 0,
          .depth = std::vector<mps::Real>(grid.cell_count(), depth),
          .momentum = std::vector<mps::Vec3>(grid.cell_count())};
}

}  // namespace

MPS_TEST_CASE("uniform resting shallow water has roundoff-scale RHS") {
  const mps::CubedSphereGrid grid(16, 6.37122e6);
  constexpr mps::Real depth = 3000.0;
  constexpr mps::Real gravity = 9.80616;
  const auto state = make_state(grid, depth);
  const auto rhs =
      mps::assemble_first_order_shallow_water_rhs(grid, state, gravity, 7.292e-5, 1.0);
  mps::Real maximum_depth_rate = 0.0;
  mps::Real maximum_momentum_rate = 0.0;
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
    maximum_depth_rate = std::max(maximum_depth_rate, std::abs(rhs.total.depth[cell]));
    maximum_momentum_rate =
        std::max(maximum_momentum_rate, mps::norm(rhs.total.momentum[cell]));
    MPS_CHECK_NEAR(mps::dot(rhs.total.momentum[cell], grid.cells()[cell].center), 0.0,
                   2.0e-13 * std::max(mps::norm(rhs.total.momentum[cell]), 1.0));
  }
  MPS_CHECK_EQ(maximum_depth_rate, 0.0);
  MPS_CHECK_NEAR(maximum_momentum_rate, 0.0, 2.0e-10);
}

MPS_TEST_CASE("unique-edge shallow-water scatter conserves global mass") {
  const mps::CubedSphereGrid grid(12, 2.0);
  auto state = make_state(grid, 3.0);
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
    state.depth[cell] += 0.2 * grid.cells()[cell].center.x;
    state.momentum[cell] =
        state.depth[cell] *
        mps::project_tangent({0.3, -0.2, 0.1}, grid.cells()[cell].center);
  }
  const auto rhs =
      mps::assemble_first_order_shallow_water_rhs(grid, state, 4.0, 0.1, 0.1);
  std::vector<mps::Real> integrated(grid.cell_count());
  mps::Real flux_scale = 0.0;
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
    integrated[cell] = rhs.total.depth[cell] * grid.cells()[cell].area_m2;
    flux_scale += std::abs(integrated[cell]);
  }
  MPS_CHECK(std::abs(mps::diagnostics::compensated_sum(integrated)) <=
            5.0e-13 * flux_scale);
}

MPS_TEST_CASE("Coriolis tendency has the documented sign and is tangent") {
  const mps::CubedSphereGrid grid(4, 1.0);
  auto state = make_state(grid, 2.0);
  const std::size_t cell = 0;
  state.momentum[cell] =
      mps::project_tangent({1.0, 2.0, 3.0}, grid.cells()[cell].center);
  const auto rhs =
      mps::assemble_first_order_shallow_water_rhs(grid, state, 3.0, 0.2, 0.1);
  const mps::Vec3 expected =
      -2.0 * mps::project_tangent(mps::cross({0.0, 0.0, 0.2}, state.momentum[cell]),
                                  grid.cells()[cell].center);
  MPS_CHECK_NEAR(mps::norm(rhs.coriolis.momentum[cell] - expected), 0.0, 1.0e-15);
  MPS_CHECK_NEAR(mps::dot(rhs.coriolis.momentum[cell], grid.cells()[cell].center), 0.0,
                 1.0e-15);
}

MPS_TEST_CASE("shallow-water RHS rejects non-finite state") {
  const mps::CubedSphereGrid grid(2, 1.0);
  auto state = make_state(grid, 2.0);
  state.depth[0] = std::numeric_limits<mps::Real>::quiet_NaN();
  MPS_CHECK_THROWS_AS(
      mps::assemble_first_order_shallow_water_rhs(grid, state, 3.0, 0.0, 0.1),
      std::runtime_error);
}

int main() { return mps::test::run_all(); }
