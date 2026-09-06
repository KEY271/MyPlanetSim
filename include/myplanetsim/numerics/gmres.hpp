#pragma once

#include <cstddef>
#include <functional>
#include <span>
#include <vector>

#include "myplanetsim/core/types.hpp"

namespace mps {

using GmresLinearOperator = std::function<void(std::span<const Real>, std::span<Real>)>;
using GmresPreconditioner = std::function<void(std::span<const Real>, std::span<Real>)>;

struct GmresOptions {
  std::size_t restart = 20;
  std::size_t maximum_iterations = 40;
  Real relative_tolerance = 1.0e-8;
  Real absolute_tolerance = 1.0e-12;
  Real breakdown_tolerance = 1.0e-14;
};

enum class GmresStatus { kConverged, kMaximumIterations, kBreakdown };

struct GmresResult {
  GmresStatus status = GmresStatus::kMaximumIterations;
  std::size_t iterations = 0;
  Real initial_residual_norm = 0.0;
  Real final_residual_norm = 0.0;
  Real relative_residual = 0.0;

  [[nodiscard]] bool converged() const noexcept {
    return status == GmresStatus::kConverged;
  }
};

struct GmresWorkspace {
  std::vector<Real> residual;
  std::vector<Real> operator_result;
  std::vector<Real> arnoldi_work;
  std::vector<Real> basis;
  std::vector<Real> preconditioned_basis;
  std::vector<Real> hessenberg;
  std::vector<Real> cosine;
  std::vector<Real> sine;
  std::vector<Real> projected_rhs;
  std::vector<Real> coefficients;
};

[[nodiscard]] GmresResult restarted_gmres(
    const GmresLinearOperator& linear_operator, std::span<const Real> right_hand_side,
    std::span<Real> solution, const GmresOptions& options, GmresWorkspace& workspace,
    const GmresPreconditioner& preconditioner = {});

}  // namespace mps
