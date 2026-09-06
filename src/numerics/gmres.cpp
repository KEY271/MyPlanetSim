#include "myplanetsim/numerics/gmres.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace mps {
namespace {

[[nodiscard]] Real dot_product(const std::span<const Real> left,
                               const std::span<const Real> right) {
  Real sum = 0.0;
  for (std::size_t index = 0; index < left.size(); ++index)
    sum += left[index] * right[index];
  return sum;
}

[[nodiscard]] Real vector_norm(const std::span<const Real> values) {
  Real result = 0.0;
  for (const Real value : values) result = std::hypot(result, value);
  return result;
}

void require_finite_vector(const std::span<const Real> values,
                           const char* const description) {
  if (!std::ranges::all_of(values,
                           [](const Real value) { return std::isfinite(value); })) {
    throw std::runtime_error(std::string(description) + " contains a non-finite value");
  }
}

void validate_inputs(const GmresLinearOperator& linear_operator,
                     const std::span<const Real> right_hand_side,
                     const std::span<const Real> solution,
                     const GmresOptions& options) {
  if (!linear_operator) throw std::invalid_argument("GMRES operator is empty");
  if (right_hand_side.empty() || solution.size() != right_hand_side.size())
    throw std::invalid_argument("GMRES vector shapes differ or are empty");
  if (options.restart == 0 || options.maximum_iterations == 0)
    throw std::invalid_argument("GMRES iteration limits must be positive");
  if (options.restart > options.maximum_iterations)
    throw std::invalid_argument("GMRES restart must not exceed maximum iterations");
  if (!(options.relative_tolerance >= 0.0) ||
      !std::isfinite(options.relative_tolerance) ||
      !(options.absolute_tolerance >= 0.0) ||
      !std::isfinite(options.absolute_tolerance) ||
      !(options.breakdown_tolerance > 0.0) ||
      !std::isfinite(options.breakdown_tolerance)) {
    throw std::invalid_argument("GMRES tolerances are invalid");
  }
  if (options.relative_tolerance == 0.0 && options.absolute_tolerance == 0.0)
    throw std::invalid_argument("GMRES requires a positive stopping tolerance");
  require_finite_vector(right_hand_side, "GMRES right-hand side");
  require_finite_vector(solution, "GMRES initial solution");
}

void resize_workspace(GmresWorkspace& workspace, const std::size_t dimension,
                      const std::size_t restart) {
  if (dimension > std::numeric_limits<std::size_t>::max() / (restart + 1))
    throw std::length_error("GMRES workspace size overflows");
  workspace.residual.resize(dimension);
  workspace.operator_result.resize(dimension);
  workspace.arnoldi_work.resize(dimension);
  workspace.basis.resize((restart + 1) * dimension);
  workspace.preconditioned_basis.resize(restart * dimension);
  workspace.hessenberg.resize((restart + 1) * restart);
  workspace.cosine.resize(restart);
  workspace.sine.resize(restart);
  workspace.projected_rhs.resize(restart + 1);
  workspace.coefficients.resize(restart);
}

[[nodiscard]] std::span<Real> vector_at(std::vector<Real>& values,
                                        const std::size_t vector_index,
                                        const std::size_t dimension) {
  return {values.data() + vector_index * dimension, dimension};
}

[[nodiscard]] std::span<const Real> vector_at(const std::vector<Real>& values,
                                              const std::size_t vector_index,
                                              const std::size_t dimension) {
  return {values.data() + vector_index * dimension, dimension};
}

[[nodiscard]] std::size_t hessenberg_offset(const std::size_t row,
                                            const std::size_t column,
                                            const std::size_t restart) {
  return row * restart + column;
}

void compute_residual(const GmresLinearOperator& linear_operator,
                      const std::span<const Real> right_hand_side,
                      const std::span<const Real> solution, GmresWorkspace& workspace) {
  linear_operator(solution, workspace.operator_result);
  require_finite_vector(workspace.operator_result, "GMRES operator result");
  for (std::size_t index = 0; index < right_hand_side.size(); ++index)
    workspace.residual[index] =
        right_hand_side[index] - workspace.operator_result[index];
}

void solve_projected_system(const std::size_t columns, const std::size_t restart,
                            GmresWorkspace& workspace) {
  std::fill(workspace.coefficients.begin(), workspace.coefficients.end(), 0.0);
  for (std::size_t reverse = 0; reverse < columns; ++reverse) {
    const std::size_t row = columns - reverse - 1;
    Real value = workspace.projected_rhs[row];
    for (std::size_t column = row + 1; column < columns; ++column) {
      value -= workspace.hessenberg[hessenberg_offset(row, column, restart)] *
               workspace.coefficients[column];
    }
    const Real diagonal = workspace.hessenberg[hessenberg_offset(row, row, restart)];
    if (diagonal == 0.0 || !std::isfinite(diagonal))
      throw std::runtime_error("GMRES projected system is singular");
    workspace.coefficients[row] = value / diagonal;
  }
}

void update_solution(const std::size_t columns, const std::size_t dimension,
                     const GmresWorkspace& workspace, const std::span<Real> solution) {
  for (std::size_t column = 0; column < columns; ++column) {
    const auto direction = vector_at(workspace.preconditioned_basis, column, dimension);
    const Real coefficient = workspace.coefficients[column];
    for (std::size_t index = 0; index < dimension; ++index)
      solution[index] += coefficient * direction[index];
  }
  require_finite_vector(solution, "GMRES solution");
}

}  // namespace

GmresResult restarted_gmres(const GmresLinearOperator& linear_operator,
                            const std::span<const Real> right_hand_side,
                            const std::span<Real> solution, const GmresOptions& options,
                            GmresWorkspace& workspace,
                            const GmresPreconditioner& preconditioner) {
  validate_inputs(linear_operator, right_hand_side, solution, options);
  const std::size_t dimension = right_hand_side.size();
  const std::size_t restart = options.restart;
  resize_workspace(workspace, dimension, restart);
  compute_residual(linear_operator, right_hand_side, solution, workspace);

  const Real right_hand_side_norm = vector_norm(right_hand_side);
  const Real relative_scale = std::max(right_hand_side_norm, Real{1});
  const Real stopping_tolerance =
      std::max(options.absolute_tolerance, options.relative_tolerance * relative_scale);
  const Real initial_residual_norm = vector_norm(workspace.residual);
  GmresResult result{.status = GmresStatus::kMaximumIterations,
                     .iterations = 0,
                     .initial_residual_norm = initial_residual_norm,
                     .final_residual_norm = initial_residual_norm,
                     .relative_residual = initial_residual_norm / relative_scale};
  if (initial_residual_norm <= stopping_tolerance) {
    result.status = GmresStatus::kConverged;
    return result;
  }

  while (result.iterations < options.maximum_iterations) {
    std::fill(workspace.hessenberg.begin(), workspace.hessenberg.end(), 0.0);
    std::fill(workspace.cosine.begin(), workspace.cosine.end(), 0.0);
    std::fill(workspace.sine.begin(), workspace.sine.end(), 0.0);
    std::fill(workspace.projected_rhs.begin(), workspace.projected_rhs.end(), 0.0);

    const Real residual_norm = vector_norm(workspace.residual);
    auto first_basis = vector_at(workspace.basis, 0, dimension);
    for (std::size_t index = 0; index < dimension; ++index)
      first_basis[index] = workspace.residual[index] / residual_norm;
    workspace.projected_rhs[0] = residual_norm;

    const std::size_t cycle_limit =
        std::min(restart, options.maximum_iterations - result.iterations);
    std::size_t completed_columns = 0;
    bool arnoldi_breakdown = false;
    for (std::size_t column = 0; column < cycle_limit; ++column) {
      const auto basis = vector_at(workspace.basis, column, dimension);
      auto direction = vector_at(workspace.preconditioned_basis, column, dimension);
      if (preconditioner) {
        preconditioner(basis, direction);
        require_finite_vector(direction, "GMRES preconditioner result");
      } else {
        std::copy(basis.begin(), basis.end(), direction.begin());
      }
      linear_operator(direction, workspace.arnoldi_work);
      require_finite_vector(workspace.arnoldi_work, "GMRES operator result");

      for (std::size_t row = 0; row <= column; ++row) {
        const auto row_basis = vector_at(workspace.basis, row, dimension);
        const Real projection = dot_product(row_basis, workspace.arnoldi_work);
        workspace.hessenberg[hessenberg_offset(row, column, restart)] = projection;
        for (std::size_t index = 0; index < dimension; ++index)
          workspace.arnoldi_work[index] -= projection * row_basis[index];
      }
      for (std::size_t row = 0; row <= column; ++row) {
        const auto row_basis = vector_at(workspace.basis, row, dimension);
        const Real correction = dot_product(row_basis, workspace.arnoldi_work);
        workspace.hessenberg[hessenberg_offset(row, column, restart)] += correction;
        for (std::size_t index = 0; index < dimension; ++index)
          workspace.arnoldi_work[index] -= correction * row_basis[index];
      }

      const Real next_norm = vector_norm(workspace.arnoldi_work);
      workspace.hessenberg[hessenberg_offset(column + 1, column, restart)] = next_norm;
      const Real column_scale = std::max(
          Real{1},
          std::abs(workspace.hessenberg[hessenberg_offset(0, column, restart)]));
      arnoldi_breakdown = next_norm <= options.breakdown_tolerance * column_scale;
      if (!arnoldi_breakdown) {
        auto next_basis = vector_at(workspace.basis, column + 1, dimension);
        for (std::size_t index = 0; index < dimension; ++index)
          next_basis[index] = workspace.arnoldi_work[index] / next_norm;
      }

      for (std::size_t row = 0; row < column; ++row) {
        const auto first = hessenberg_offset(row, column, restart);
        const auto second = hessenberg_offset(row + 1, column, restart);
        const Real upper = workspace.hessenberg[first];
        const Real lower = workspace.hessenberg[second];
        workspace.hessenberg[first] =
            workspace.cosine[row] * upper + workspace.sine[row] * lower;
        workspace.hessenberg[second] =
            -workspace.sine[row] * upper + workspace.cosine[row] * lower;
      }

      const auto diagonal = hessenberg_offset(column, column, restart);
      const auto subdiagonal = hessenberg_offset(column + 1, column, restart);
      const Real rotation_norm =
          std::hypot(workspace.hessenberg[diagonal], workspace.hessenberg[subdiagonal]);
      if (rotation_norm == 0.0 || !std::isfinite(rotation_norm)) {
        result.status = GmresStatus::kBreakdown;
        return result;
      }
      workspace.cosine[column] = workspace.hessenberg[diagonal] / rotation_norm;
      workspace.sine[column] = workspace.hessenberg[subdiagonal] / rotation_norm;
      workspace.hessenberg[diagonal] = rotation_norm;
      workspace.hessenberg[subdiagonal] = 0.0;

      const Real projected = workspace.projected_rhs[column];
      workspace.projected_rhs[column] = workspace.cosine[column] * projected;
      workspace.projected_rhs[column + 1] = -workspace.sine[column] * projected;
      ++result.iterations;
      completed_columns = column + 1;

      if (std::abs(workspace.projected_rhs[column + 1]) <= stopping_tolerance ||
          arnoldi_breakdown)
        break;
    }

    solve_projected_system(completed_columns, restart, workspace);
    update_solution(completed_columns, dimension, workspace, solution);
    compute_residual(linear_operator, right_hand_side, solution, workspace);
    result.final_residual_norm = vector_norm(workspace.residual);
    result.relative_residual = result.final_residual_norm / relative_scale;
    if (result.final_residual_norm <= stopping_tolerance) {
      result.status = GmresStatus::kConverged;
      return result;
    }
    if (arnoldi_breakdown) {
      result.status = GmresStatus::kBreakdown;
      return result;
    }
  }

  result.status = GmresStatus::kMaximumIterations;
  return result;
}

}  // namespace mps
