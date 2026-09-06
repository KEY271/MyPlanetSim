#include "myplanetsim/numerics/helmholtz_solver.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace mps {
namespace {

void validate_field(const CubedSphereGrid& grid, const std::span<const Real> values,
                    const char* const description) {
  if (values.size() != grid.cell_count())
    throw std::invalid_argument(std::string(description) +
                                " shape does not match grid");
  if (!std::ranges::all_of(values,
                           [](const Real value) { return std::isfinite(value); })) {
    throw std::invalid_argument(std::string(description) +
                                " contains a non-finite value");
  }
}

}  // namespace

FiniteVolumeHelmholtzOperator::FiniteVolumeHelmholtzOperator(
    const CubedSphereGrid& grid, const Real coefficient_m2)
    : grid_(&grid),
      coefficient_m2_(coefficient_m2),
      inverse_diagonal_(grid.cell_count(), 1.0) {
  if (!(coefficient_m2 >= 0.0) || !std::isfinite(coefficient_m2))
    throw std::invalid_argument(
        "Helmholtz coefficient must be finite and non-negative");
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
    Real laplacian_diagonal = 0.0;
    for (const auto& edge : grid.cell_cache()[cell].edges) {
      laplacian_diagonal +=
          grid.edges()[edge.edge].length_m /
          (grid.cells()[cell].area_m2 * grid.edge_cache()[edge.edge].center_distance_m);
    }
    inverse_diagonal_[cell] = 1.0 / (1.0 + coefficient_m2 * laplacian_diagonal);
  }
}

void FiniteVolumeHelmholtzOperator::apply_laplacian(
    const std::span<const Real> input, const std::span<Real> output) const {
  validate_field(*grid_, input, "Helmholtz input");
  if (output.size() != grid_->cell_count())
    throw std::invalid_argument("Helmholtz output shape does not match grid");
  if (input.data() == output.data())
    throw std::invalid_argument("Helmholtz input and output must not alias");
  std::fill(output.begin(), output.end(), 0.0);
  for (const auto& edge : grid_->edges()) {
    const auto& cached = grid_->edge_cache()[edge.id];
    const Real flux = edge.length_m *
                      (input[cached.right_cell] - input[cached.left_cell]) /
                      cached.center_distance_m;
    output[cached.left_cell] += flux / grid_->cells()[cached.left_cell].area_m2;
    output[cached.right_cell] -= flux / grid_->cells()[cached.right_cell].area_m2;
  }
}

void FiniteVolumeHelmholtzOperator::apply(const std::span<const Real> input,
                                          const std::span<Real> output) const {
  apply_laplacian(input, output);
  for (std::size_t cell = 0; cell < input.size(); ++cell)
    output[cell] = input[cell] - coefficient_m2_ * output[cell];
}

void FiniteVolumeHelmholtzOperator::apply_jacobi_preconditioner(
    const std::span<const Real> input, const std::span<Real> output) const {
  validate_field(*grid_, input, "Helmholtz preconditioner input");
  if (output.size() != grid_->cell_count())
    throw std::invalid_argument(
        "Helmholtz preconditioner output shape does not match grid");
  for (std::size_t cell = 0; cell < input.size(); ++cell)
    output[cell] = inverse_diagonal_[cell] * input[cell];
}

GmresResult solve_finite_volume_helmholtz(
    const FiniteVolumeHelmholtzOperator& helmholtz,
    const std::span<const Real> right_hand_side, const std::span<Real> solution,
    const GmresOptions& options, GmresWorkspace& workspace) {
  const GmresLinearOperator linear_operator =
      [&helmholtz](const std::span<const Real> input, const std::span<Real> output) {
        helmholtz.apply(input, output);
      };
  const GmresPreconditioner preconditioner =
      [&helmholtz](const std::span<const Real> input, const std::span<Real> output) {
        helmholtz.apply_jacobi_preconditioner(input, output);
      };
  return restarted_gmres(linear_operator, right_hand_side, solution, options, workspace,
                         preconditioner);
}

}  // namespace mps
