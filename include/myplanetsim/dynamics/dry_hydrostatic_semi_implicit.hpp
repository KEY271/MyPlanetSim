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
};

struct DryHydrostaticExternalSolveResult {
  GmresResult linear;
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

}  // namespace mps
