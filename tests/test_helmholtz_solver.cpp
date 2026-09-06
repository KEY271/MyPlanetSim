#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <vector>

#include "myplanetsim/numerics/helmholtz_solver.hpp"
#include "support/test.hpp"

MPS_TEST_CASE("finite-volume Helmholtz preserves constants") {
  const mps::CubedSphereGrid grid(4, 2.0);
  const mps::FiniteVolumeHelmholtzOperator helmholtz(grid, 3.0);
  const std::vector<mps::Real> constant(grid.cell_count(), 7.0);
  std::vector<mps::Real> result(grid.cell_count());
  helmholtz.apply(constant, result);
  for (const auto value : result) MPS_CHECK_NEAR(value, 7.0, 1.0e-14);
}

MPS_TEST_CASE("finite-volume Laplacian has zero area-weighted integral") {
  const mps::CubedSphereGrid grid(6, 3.0);
  const mps::FiniteVolumeHelmholtzOperator helmholtz(grid, 1.0);
  std::vector<mps::Real> field(grid.cell_count());
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell)
    field[cell] = grid.cells()[cell].center.x +
                  0.3 * grid.cells()[cell].center.y * grid.cells()[cell].center.z;
  std::vector<mps::Real> laplacian(grid.cell_count());
  helmholtz.apply_laplacian(field, laplacian);
  mps::Real integral = 0.0;
  mps::Real absolute_integral = 0.0;
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
    integral += grid.cells()[cell].area_m2 * laplacian[cell];
    absolute_integral += std::abs(grid.cells()[cell].area_m2 * laplacian[cell]);
  }
  MPS_CHECK_NEAR(integral, 0.0, 2.0e-14 * absolute_integral);
}

MPS_TEST_CASE("GMRES solves the finite-volume Helmholtz system") {
  const mps::CubedSphereGrid grid(5, 1.0);
  const mps::FiniteVolumeHelmholtzOperator helmholtz(grid, 0.4);
  std::vector<mps::Real> expected(grid.cell_count());
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell)
    expected[cell] =
        1.0 + 0.2 * grid.cells()[cell].center.x - 0.1 * grid.cells()[cell].center.z;
  std::vector<mps::Real> right_hand_side(grid.cell_count());
  helmholtz.apply(expected, right_hand_side);
  std::vector<mps::Real> solution(grid.cell_count(), 0.0);
  mps::GmresWorkspace workspace;
  const auto result =
      mps::solve_finite_volume_helmholtz(helmholtz, right_hand_side, solution,
                                         {.restart = 20,
                                          .maximum_iterations = 80,
                                          .relative_tolerance = 1.0e-11,
                                          .absolute_tolerance = 1.0e-13},
                                         workspace);
  MPS_CHECK(result.converged());
  MPS_CHECK(result.iterations < 80);
  mps::Real maximum_error = 0.0;
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell)
    maximum_error = std::max(maximum_error, std::abs(solution[cell] - expected[cell]));
  MPS_CHECK(maximum_error < 1.0e-9);
}

MPS_TEST_CASE("Helmholtz rejects invalid coefficient and shape") {
  const mps::CubedSphereGrid grid(2, 1.0);
  MPS_CHECK_THROWS_AS(mps::FiniteVolumeHelmholtzOperator(grid, -1.0),
                      std::invalid_argument);
  const mps::FiniteVolumeHelmholtzOperator helmholtz(grid, 1.0);
  std::vector<mps::Real> input(grid.cell_count());
  std::vector<mps::Real> short_output(grid.cell_count() - 1);
  MPS_CHECK_THROWS_AS(helmholtz.apply(input, short_output), std::invalid_argument);
  MPS_CHECK_THROWS_AS(helmholtz.apply(input, input), std::invalid_argument);
}

int main() { return mps::test::run_all(); }
