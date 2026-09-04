#include <algorithm>
#include <cmath>
#include <vector>

#include "myplanetsim/numerics/spherical_operators.hpp"
#include "support/test.hpp"

namespace {

[[nodiscard]] mps::Vec3 direct_scalar_gradient(const mps::CubedSphereGrid& grid,
                                               const std::vector<double>& values,
                                               const std::size_t cell) {
  const auto& cached = grid.cell_cache()[cell];
  double a00 = 0.0;
  double a01 = 0.0;
  double a11 = 0.0;
  double b0 = 0.0;
  double b1 = 0.0;
  for (const auto& edge : cached.edges) {
    const double x = edge.neighbor_coordinates_m.alpha;
    const double y = edge.neighbor_coordinates_m.beta;
    const double weight = 1.0 / std::max(x * x + y * y, 1.0e-300);
    const double difference = values[edge.neighbor] - values[cell];
    a00 += weight * x * x;
    a01 += weight * x * y;
    a11 += weight * y * y;
    b0 += weight * x * difference;
    b1 += weight * y * difference;
  }
  const double determinant = a00 * a11 - a01 * a01;
  const double alpha = (a11 * b0 - a01 * b1) / determinant;
  const double beta = (a00 * b1 - a01 * b0) / determinant;
  return alpha * cached.basis.alpha + beta * cached.basis.beta;
}

}  // namespace

MPS_TEST_CASE("unique edge scatter cancels globally") {
  const mps::CubedSphereGrid grid(8, 2.0);
  std::vector<double> flux(grid.edge_count());
  for (std::size_t edge = 0; edge < flux.size(); ++edge) {
    flux[edge] = std::sin(static_cast<double>(edge));
  }
  const auto divergence = mps::finite_volume_divergence(grid, flux);
  double integral = 0.0;
  for (std::size_t cell = 0; cell < divergence.size(); ++cell) {
    integral += divergence[cell] * grid.cells()[cell].area_m2;
  }
  MPS_CHECK_NEAR(integral, 0.0, 2.0e-13 * static_cast<double>(flux.size()));
}

MPS_TEST_CASE("least squares gradient annihilates constants") {
  const mps::CubedSphereGrid grid(12, 3.0);
  std::vector<double> constant(grid.cell_count(), 4.25);
  const auto gradient = mps::least_squares_gradient(grid, constant);
  for (const auto value : gradient) {
    MPS_CHECK_EQ(mps::norm_squared(value), 0.0);
  }
  const auto laplacian = mps::finite_volume_laplacian(grid, constant);
  for (const auto value : laplacian) {
    MPS_CHECK_EQ(value, 0.0);
  }
}

MPS_TEST_CASE("manufactured gradient error decreases with resolution") {
  double previous = 1.0e100;
  for (const mps::Index n : {8, 16, 32}) {
    const mps::CubedSphereGrid grid(n, 1.0);
    std::vector<double> field(grid.cell_count());
    for (std::size_t i = 0; i < field.size(); ++i) {
      field[i] = grid.cells()[i].center.x;
    }
    const auto gradient = mps::least_squares_gradient(grid, field);
    double error = 0.0;
    for (std::size_t i = 0; i < field.size(); ++i) {
      const auto exact = mps::project_tangent({1.0, 0.0, 0.0}, grid.cells()[i].center);
      error += mps::norm_squared(gradient[i] - exact) * grid.cells()[i].area_m2;
    }
    error = std::sqrt(error / grid.total_area_m2());
    MPS_CHECK(error < previous);
    previous = error;
  }
}

MPS_TEST_CASE("precomputed least squares weights match the direct solve") {
  const mps::CubedSphereGrid grid(9, 3.0);
  std::vector<double> field(grid.cell_count());
  for (std::size_t cell = 0; cell < field.size(); ++cell) {
    const auto center = grid.cells()[cell].center;
    field[cell] = std::sin(1.3 * center.x) + 0.25 * center.y * center.z;
  }

  const auto gradient = mps::least_squares_gradient(grid, field);
  double maximum_error = 0.0;
  double maximum_reference = 0.0;
  for (std::size_t cell = 0; cell < field.size(); ++cell) {
    const auto reference = direct_scalar_gradient(grid, field, cell);
    maximum_error = std::max(maximum_error, mps::norm(gradient[cell] - reference));
    maximum_reference = std::max(maximum_reference, mps::norm(reference));
  }
  MPS_CHECK(maximum_error <= 1.0e-12 * maximum_reference);
}

int main() { return mps::test::run_all(); }
