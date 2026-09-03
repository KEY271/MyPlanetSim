#include <algorithm>
#include <cmath>
#include <vector>

#include "myplanetsim/numerics/spherical_operators.hpp"
#include "support/test.hpp"

namespace {

[[nodiscard]] std::vector<mps::Vec3> solid_rotation(const mps::CubedSphereGrid& grid,
                                                    const mps::Vec3 omega) {
  std::vector<mps::Vec3> velocity(grid.cell_count());
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
    velocity[cell] = grid.radius_m() * mps::cross(omega, grid.cells()[cell].center);
  }
  return velocity;
}

[[nodiscard]] mps::Real reconstruction_error(const mps::Index n) {
  const mps::CubedSphereGrid grid(n, 2.0);
  constexpr mps::Vec3 omega{0.2, -0.1, 0.3};
  const auto velocity = solid_rotation(grid, omega);
  const auto gradients = mps::least_squares_vector_gradient(grid, velocity);
  mps::Real squared_error = 0.0;
  for (const auto& edge : grid.edges()) {
    const std::size_t left = grid.cell_index(edge.left_cell);
    const mps::Vec3 reconstructed = mps::reconstruct_tangent_vector(
        grid, left, velocity[left], edge.center, gradients[left]);
    const mps::Vec3 exact = grid.radius_m() * mps::cross(omega, edge.center);
    squared_error += mps::norm_squared(reconstructed - exact);
  }
  return std::sqrt(squared_error / static_cast<mps::Real>(grid.edge_count()));
}

}  // namespace

MPS_TEST_CASE("edge tangent basis is orthonormal and consistently oriented") {
  const mps::CubedSphereGrid grid(8, 1.0);
  for (const auto& edge : grid.edges()) {
    const auto basis = mps::edge_tangent_basis(edge);
    MPS_CHECK_NEAR(mps::norm(basis.normal), 1.0, 5.0e-15);
    MPS_CHECK_NEAR(mps::norm(basis.tangent), 1.0, 5.0e-15);
    MPS_CHECK_NEAR(mps::dot(basis.normal, basis.tangent), 0.0, 5.0e-15);
    MPS_CHECK(mps::dot(basis.normal, edge.outward_normal_from_left) > 0.0);
  }
}

MPS_TEST_CASE("tangent-vector reconstruction converges for solid rotation") {
  const mps::Real coarse = reconstruction_error(8);
  const mps::Real fine = reconstruction_error(16);
  MPS_CHECK(fine < 0.35 * coarse);
}

MPS_TEST_CASE("shallow-water vector operators preserve zero fields") {
  const mps::CubedSphereGrid grid(8, 3.0);
  const std::vector<mps::Vec3> zero(grid.cell_count());
  const auto divergence = mps::finite_volume_vector_divergence(grid, zero);
  const auto curl = mps::finite_volume_curl(grid, zero);
  const auto laplacian = mps::finite_volume_vector_laplacian(grid, zero);
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
    MPS_CHECK_EQ(divergence[cell], 0.0);
    MPS_CHECK_EQ(curl[cell], 0.0);
    MPS_CHECK_EQ(mps::norm(laplacian[cell]), 0.0);
  }
}

MPS_TEST_CASE("potential vorticity includes relative and planetary vorticity") {
  const mps::CubedSphereGrid grid(24, 2.0);
  constexpr mps::Vec3 omega{0.0, 0.0, 0.25};
  const auto velocity = solid_rotation(grid, omega);
  const std::vector<mps::Real> depth(grid.cell_count(), 4.0);
  const auto pv = mps::shallow_water_potential_vorticity(grid, depth, velocity, omega);
  mps::Real maximum_error = 0.0;
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
    const mps::Real exact =
        4.0 * mps::dot(omega, grid.cells()[cell].center) / depth[cell];
    maximum_error = std::max(maximum_error, std::abs(pv[cell] - exact));
  }
  MPS_CHECK(maximum_error < 0.015);
}

int main() { return mps::test::run_all(); }
