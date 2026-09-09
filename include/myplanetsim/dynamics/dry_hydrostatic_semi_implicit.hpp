#pragma once

#include <span>
#include <vector>

#include "myplanetsim/dynamics/dry_hydrostatic_fast_operator.hpp"
#include "myplanetsim/numerics/gmres.hpp"

namespace mps {

struct DryHydrostaticExternalModeOperator {
  Real phase_speed_m_s = 0.0;
  std::vector<Real> momentum_weights;
};

struct DryHydrostaticSemiImplicitWorkspace {
  DryHydrostaticFastOperatorWorkspace fast_operator;
  GmresWorkspace gmres;
  std::vector<Vec3> column_momentum;
  std::vector<Real> column_divergence;
  std::vector<Real> helmholtz_right_hand_side;
  std::vector<Real> surface_pressure_correction;
  std::vector<Vec3> surface_pressure_gradient;
  std::vector<Vec3> helmholtz_gradient;
  std::vector<Real> helmholtz_divergence;
  DryHydrostaticFastPerturbation correction;
  DryHydrostaticFastTendency external_tendency;
  DryHydrostaticFastPerturbation scalar_right_hand_side;
  DryHydrostaticFastTendency scalar_force;
  DryHydrostaticFastPerturbation momentum_perturbation;
  DryHydrostaticFastTendency momentum_scalar_tendency;
  std::vector<Vec3> effective_momentum;
  std::vector<Vec3> modal_momentum;
  std::vector<Vec3> modal_momentum_right_hand_side;
  std::vector<Vec3> modal_delta;
  DryHydrostaticModalBatchWorkspace modal_batch;
  std::vector<Real> modal_divergence;
  std::vector<Real> modal_solution;
  std::vector<Real> vector_divergence_edge_flux;
  std::vector<Real> helmholtz_laplacian_diagonal;
  std::vector<Real> helmholtz_inverse_diagonal;
  std::vector<unsigned char> selected_mode_mask;
  const CubedSphereGrid* active_modal_grid = nullptr;
  Real active_modal_coefficient_m2 = 0.0;
  GmresLinearOperator modal_helmholtz_operator;
  GmresPreconditioner modal_helmholtz_preconditioner;
};

struct DryHydrostaticExternalSolveResult {
  GmresResult linear;
  Real equation_residual_norm = 0.0;
};

struct DryHydrostaticModalSolveResult {
  bool all_converged = false;
  std::size_t selected_modes = 0;
  std::size_t linear_iterations_total = 0;
  std::size_t linear_iterations_maximum = 0;
  Real linear_relative_residual_maximum = 0.0;
  Real equation_residual_norm = 0.0;
};

[[nodiscard]] DryHydrostaticExternalModeOperator
make_dry_hydrostatic_external_mode_operator(
    const DryHydrostaticReferenceColumn& reference,
    const DryHydrostaticVerticalModes& modes);

void apply_dry_hydrostatic_external_mode_operator(
    const CubedSphereGrid& grid, const PlanetParameters& planet,
    const DryHydrostaticFastOperator& fast_operator,
    const DryHydrostaticExternalModeOperator& external_operator,
    const DryHydrostaticFastPerturbation& perturbation,
    DryHydrostaticFastTendency& result, DryHydrostaticFastOperatorWorkspace& workspace);

// Solve [I - implicit_time_s * L_external] correction = right_hand_side.
// The scalar Schur complement is a two-dimensional cell-centred Helmholtz
// problem using the same centred divergence and least-squares gradient as L.
[[nodiscard]] DryHydrostaticExternalSolveResult
solve_dry_hydrostatic_external_mode_correction(
    const CubedSphereGrid& grid, const PlanetParameters& planet,
    const DryHydrostaticFastOperator& fast_operator,
    const DryHydrostaticExternalModeOperator& external_operator, Real implicit_time_s,
    const DryHydrostaticFastPerturbation& right_hand_side, const GmresOptions& options,
    DryHydrostaticFastPerturbation& correction,
    DryHydrostaticSemiImplicitWorkspace& workspace);

// Apply the exact vertical modal factorization to the selected wave modes. Modes
// below the selection threshold retain the identity correction, leaving their
// small-Courant response to the outer nonlinear iteration.
[[nodiscard]] DryHydrostaticModalSolveResult solve_dry_hydrostatic_modal_correction(
    const CubedSphereGrid& grid, const PlanetParameters& planet,
    const DryHydrostaticFastOperator& fast_operator,
    const DryHydrostaticVerticalModes& modes,
    std::span<const std::size_t> selected_modes, Real implicit_time_s,
    const DryHydrostaticFastPerturbation& right_hand_side, const GmresOptions& options,
    DryHydrostaticFastPerturbation& correction,
    DryHydrostaticSemiImplicitWorkspace& workspace);

}  // namespace mps
