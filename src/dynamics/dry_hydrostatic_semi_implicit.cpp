#include "myplanetsim/dynamics/dry_hydrostatic_semi_implicit.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "myplanetsim/numerics/helmholtz_solver.hpp"
#include "myplanetsim/numerics/spherical_operators.hpp"

namespace mps {
namespace {

[[nodiscard]] Real scaled_equation_residual(
    const DryHydrostaticFastPerturbation& actual,
    const DryHydrostaticFastPerturbation& expected) {
  Real residual = 0.0;
  for (std::size_t cell = 0; cell < actual.surface_pressure_pa.size(); ++cell) {
    residual = std::max(
        residual, std::abs(actual.surface_pressure_pa[cell] -
                           expected.surface_pressure_pa[cell]) /
                      std::max(1.0, std::abs(expected.surface_pressure_pa[cell])));
  }
  for (std::size_t n = 0; n < actual.horizontal_momentum_mass_kg_m_s.size(); ++n) {
    residual = std::max(
        residual, norm(actual.horizontal_momentum_mass_kg_m_s[n] -
                       expected.horizontal_momentum_mass_kg_m_s[n]) /
                      std::max(1.0, norm(expected.horizontal_momentum_mass_kg_m_s[n])));
    residual = std::max(
        residual,
        std::abs(actual.potential_temperature_mass_k_kg_m2[n] -
                 expected.potential_temperature_mass_k_kg_m2[n]) /
            std::max(1.0, std::abs(expected.potential_temperature_mass_k_kg_m2[n])));
  }
  return residual;
}

}  // namespace

DryHydrostaticExternalModeOperator make_dry_hydrostatic_external_mode_operator(
    const DryHydrostaticReferenceColumn& reference,
    const DryHydrostaticVerticalModes& modes) {
  const auto levels = reference.geometry.air_mass_kg_m2.size();
  if (levels == 0 || modes.mode_count() == 0 || modes.levels != levels ||
      !(modes.phase_speed_m_s.front() > 0.0) ||
      !std::isfinite(modes.phase_speed_m_s.front()))
    throw std::invalid_argument("external dry mode is unavailable");
  Real column_mass = 0.0;
  for (const auto mass : reference.geometry.air_mass_kg_m2) column_mass += mass;
  if (!(column_mass > 0.0) || !std::isfinite(column_mass))
    throw std::invalid_argument("reference column mass is invalid");
  DryHydrostaticExternalModeOperator result{
      .phase_speed_m_s = modes.phase_speed_m_s.front(),
      .momentum_weights = std::vector<Real>(levels)};
  for (std::size_t level = 0; level < levels; ++level)
    result.momentum_weights[level] =
        reference.geometry.air_mass_kg_m2[level] / column_mass;
  return result;
}

void apply_dry_hydrostatic_external_mode_operator(
    const CubedSphereGrid& grid, const PlanetParameters& planet,
    const DryHydrostaticFastOperator& fast_operator,
    const DryHydrostaticExternalModeOperator& external_operator,
    const DryHydrostaticFastPerturbation& perturbation,
    DryHydrostaticFastTendency& result,
    DryHydrostaticFastOperatorWorkspace& workspace) {
  if (external_operator.momentum_weights.size() != fast_operator.levels)
    throw std::invalid_argument("external mode and fast operator shapes differ");
  apply_dry_hydrostatic_fast_operator(grid, planet, fast_operator, perturbation, result,
                                      workspace);
  workspace.pressure_gradient.resize(grid.cell_count());
  least_squares_gradient(grid, perturbation.surface_pressure_pa,
                         workspace.pressure_gradient);
  const Real acceleration_coefficient = external_operator.phase_speed_m_s *
                                        external_operator.phase_speed_m_s /
                                        planet.gravity_m_s2;
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
    for (std::size_t level = 0; level < fast_operator.levels; ++level) {
      const auto n = dry_hydrostatic_offset(cell, level, fast_operator.levels);
      result.tendency.momentum[n] = -external_operator.momentum_weights[level] *
                                    acceleration_coefficient *
                                    workspace.pressure_gradient[cell];
    }
  }
}

DryHydrostaticExternalSolveResult solve_dry_hydrostatic_external_mode_correction(
    const CubedSphereGrid& grid, const PlanetParameters& planet,
    const DryHydrostaticFastOperator& fast_operator,
    const DryHydrostaticExternalModeOperator& external_operator,
    const Real implicit_time_s, const DryHydrostaticFastPerturbation& right_hand_side,
    const GmresOptions& options, DryHydrostaticFastPerturbation& correction,
    DryHydrostaticSemiImplicitWorkspace& workspace) {
  if (!(implicit_time_s > 0.0) || !std::isfinite(implicit_time_s))
    throw std::invalid_argument("implicit dry correction time must be positive");
  const auto cells = grid.cell_count();
  const auto levels = fast_operator.levels;
  const auto volume = cells * levels;
  if (right_hand_side.surface_pressure_pa.size() != cells ||
      right_hand_side.horizontal_momentum_mass_kg_m_s.size() != volume ||
      right_hand_side.potential_temperature_mass_k_kg_m2.size() != volume ||
      external_operator.momentum_weights.size() != levels)
    throw std::invalid_argument("external correction right-hand side shape mismatch");

  workspace.column_momentum.assign(cells, {});
  for (std::size_t cell = 0; cell < cells; ++cell) {
    for (std::size_t level = 0; level < levels; ++level) {
      workspace.column_momentum[cell] =
          workspace.column_momentum[cell] +
          right_hand_side.horizontal_momentum_mass_kg_m_s[dry_hydrostatic_offset(
              cell, level, levels)];
    }
  }
  workspace.column_divergence =
      finite_volume_vector_divergence(grid, workspace.column_momentum);
  workspace.helmholtz_right_hand_side.resize(cells);
  for (std::size_t cell = 0; cell < cells; ++cell) {
    workspace.helmholtz_right_hand_side[cell] =
        right_hand_side.surface_pressure_pa[cell] -
        implicit_time_s * planet.gravity_m_s2 * workspace.column_divergence[cell];
  }

  const Real coefficient_m2 = implicit_time_s * implicit_time_s *
                              external_operator.phase_speed_m_s *
                              external_operator.phase_speed_m_s;
  const FiniteVolumeHelmholtzOperator preconditioner(grid, coefficient_m2);
  workspace.helmholtz_gradient.resize(cells);
  workspace.helmholtz_divergence.resize(cells);
  const GmresLinearOperator helmholtz = [&](const std::span<const Real> input,
                                            const std::span<Real> output) {
    least_squares_gradient(grid, input, workspace.helmholtz_gradient);
    const auto divergence =
        finite_volume_vector_divergence(grid, workspace.helmholtz_gradient);
    std::copy(divergence.begin(), divergence.end(),
              workspace.helmholtz_divergence.begin());
    for (std::size_t cell = 0; cell < cells; ++cell)
      output[cell] = input[cell] - coefficient_m2 * divergence[cell];
  };
  const GmresPreconditioner jacobi = [&](const std::span<const Real> input,
                                         const std::span<Real> output) {
    preconditioner.apply_jacobi_preconditioner(input, output);
  };
  workspace.surface_pressure_correction = right_hand_side.surface_pressure_pa;
  const auto linear = restarted_gmres(helmholtz, workspace.helmholtz_right_hand_side,
                                      workspace.surface_pressure_correction, options,
                                      workspace.gmres, jacobi);
  if (!linear.converged()) return {.linear = linear};

  correction.surface_pressure_pa = workspace.surface_pressure_correction;
  correction.horizontal_momentum_mass_kg_m_s =
      right_hand_side.horizontal_momentum_mass_kg_m_s;
  correction.potential_temperature_mass_k_kg_m2 =
      right_hand_side.potential_temperature_mass_k_kg_m2;
  workspace.surface_pressure_gradient.resize(cells);
  least_squares_gradient(grid, correction.surface_pressure_pa,
                         workspace.surface_pressure_gradient);
  const Real acceleration_coefficient = external_operator.phase_speed_m_s *
                                        external_operator.phase_speed_m_s /
                                        planet.gravity_m_s2;
  for (std::size_t cell = 0; cell < cells; ++cell) {
    for (std::size_t level = 0; level < levels; ++level) {
      const auto n = dry_hydrostatic_offset(cell, level, levels);
      correction.horizontal_momentum_mass_kg_m_s[n] =
          correction.horizontal_momentum_mass_kg_m_s[n] -
          implicit_time_s * external_operator.momentum_weights[level] *
              acceleration_coefficient * workspace.surface_pressure_gradient[cell];
    }
  }

  DryHydrostaticFastPerturbation pressure_momentum_correction = correction;
  std::fill(pressure_momentum_correction.potential_temperature_mass_k_kg_m2.begin(),
            pressure_momentum_correction.potential_temperature_mass_k_kg_m2.end(), 0.0);
  apply_dry_hydrostatic_external_mode_operator(
      grid, planet, fast_operator, external_operator, pressure_momentum_correction,
      workspace.external_tendency, workspace.fast_operator);
  for (std::size_t n = 0; n < volume; ++n) {
    correction.potential_temperature_mass_k_kg_m2[n] +=
        implicit_time_s *
        workspace.external_tendency.tendency.potential_temperature_mass[n];
  }

  apply_dry_hydrostatic_external_mode_operator(
      grid, planet, fast_operator, external_operator, correction,
      workspace.external_tendency, workspace.fast_operator);
  DryHydrostaticFastPerturbation equation = correction;
  for (std::size_t cell = 0; cell < cells; ++cell) {
    equation.surface_pressure_pa[cell] -=
        implicit_time_s * workspace.external_tendency.surface_pressure_pa_s[cell];
  }
  for (std::size_t n = 0; n < volume; ++n) {
    equation.horizontal_momentum_mass_kg_m_s[n] =
        equation.horizontal_momentum_mass_kg_m_s[n] -
        implicit_time_s * workspace.external_tendency.tendency.momentum[n];
    equation.potential_temperature_mass_k_kg_m2[n] -=
        implicit_time_s *
        workspace.external_tendency.tendency.potential_temperature_mass[n];
  }
  return {
      .linear = linear,
      .equation_residual_norm = scaled_equation_residual(equation, right_hand_side)};
}

}  // namespace mps
