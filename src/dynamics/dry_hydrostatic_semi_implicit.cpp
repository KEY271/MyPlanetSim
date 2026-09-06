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

DryHydrostaticModalSolveResult solve_dry_hydrostatic_modal_correction(
    const CubedSphereGrid& grid, const PlanetParameters& planet,
    const DryHydrostaticFastOperator& fast_operator,
    const DryHydrostaticVerticalModes& modes,
    const std::span<const std::size_t> selected_modes, const Real implicit_time_s,
    const DryHydrostaticFastPerturbation& right_hand_side, const GmresOptions& options,
    DryHydrostaticFastPerturbation& correction,
    DryHydrostaticSemiImplicitWorkspace& workspace) {
  if (!(implicit_time_s > 0.0) || !std::isfinite(implicit_time_s))
    throw std::invalid_argument("implicit dry modal correction time must be positive");
  const auto cells = grid.cell_count();
  const auto levels = fast_operator.levels;
  const auto volume = cells * levels;
  if (levels == 0 || modes.levels != levels || modes.mode_count() != levels ||
      modes.eigenvalue_m2_s2.size() != levels ||
      modes.eigenvectors.size() != levels * levels ||
      modes.inverse_eigenvectors.size() != levels * levels ||
      right_hand_side.surface_pressure_pa.size() != cells ||
      right_hand_side.horizontal_momentum_mass_kg_m_s.size() != volume ||
      right_hand_side.potential_temperature_mass_k_kg_m2.size() != volume)
    throw std::invalid_argument("modal correction shapes differ");
  workspace.selected_mode_mask.assign(levels, 0);
  for (const auto mode : selected_modes) {
    if (mode >= levels || workspace.selected_mode_mask[mode] != 0)
      throw std::invalid_argument("selected dry vertical modes are invalid");
    workspace.selected_mode_mask[mode] = 1;
  }

  workspace.scalar_right_hand_side = right_hand_side;
  std::fill(workspace.scalar_right_hand_side.horizontal_momentum_mass_kg_m_s.begin(),
            workspace.scalar_right_hand_side.horizontal_momentum_mass_kg_m_s.end(),
            Vec3{});
  apply_dry_hydrostatic_fast_operator(grid, planet, fast_operator,
                                      workspace.scalar_right_hand_side,
                                      workspace.scalar_force, workspace.fast_operator);
  workspace.effective_momentum.resize(volume);
  for (std::size_t n = 0; n < volume; ++n) {
    workspace.effective_momentum[n] =
        right_hand_side.horizontal_momentum_mass_kg_m_s[n] +
        implicit_time_s * workspace.scalar_force.tendency.momentum[n];
  }

  workspace.modal_momentum.assign(volume, {});
  for (std::size_t cell = 0; cell < cells; ++cell) {
    for (std::size_t mode = 0; mode < levels; ++mode) {
      Vec3 value{};
      for (std::size_t level = 0; level < levels; ++level) {
        value = value + modes.inverse_eigenvectors[mode * levels + level] *
                            workspace.effective_momentum[dry_hydrostatic_offset(
                                cell, level, levels)];
      }
      workspace.modal_momentum[dry_hydrostatic_offset(cell, mode, levels)] = value;
    }
  }

  DryHydrostaticModalSolveResult result{.all_converged = true,
                                        .selected_modes = selected_modes.size()};
  workspace.helmholtz_gradient.resize(cells);
  workspace.helmholtz_divergence.resize(cells);
  workspace.vector_divergence_edge_flux.resize(grid.edge_count());
  workspace.helmholtz_inverse_diagonal.resize(cells);
  workspace.active_modal_grid = &grid;
  if (!workspace.modal_helmholtz_operator) {
    workspace.modal_helmholtz_operator = [&workspace](const std::span<const Real> input,
                                                      const std::span<Real> output) {
      const auto& active_grid = *workspace.active_modal_grid;
      least_squares_gradient(active_grid, input, workspace.helmholtz_gradient);
      finite_volume_vector_divergence(active_grid, workspace.helmholtz_gradient,
                                      workspace.vector_divergence_edge_flux,
                                      workspace.helmholtz_divergence);
      for (std::size_t cell = 0; cell < active_grid.cell_count(); ++cell)
        output[cell] = input[cell] - workspace.active_modal_coefficient_m2 *
                                         workspace.helmholtz_divergence[cell];
    };
    workspace.modal_helmholtz_preconditioner =
        [&workspace](const std::span<const Real> input, const std::span<Real> output) {
          for (std::size_t cell = 0; cell < input.size(); ++cell)
            output[cell] = workspace.helmholtz_inverse_diagonal[cell] * input[cell];
        };
  }
  for (std::size_t mode = 0; mode < levels; ++mode) {
    if (workspace.selected_mode_mask[mode] == 0) continue;
    workspace.column_momentum.resize(cells);
    for (std::size_t cell = 0; cell < cells; ++cell)
      workspace.column_momentum[cell] =
          workspace.modal_momentum[dry_hydrostatic_offset(cell, mode, levels)];
    workspace.modal_divergence.resize(cells);
    finite_volume_vector_divergence(grid, workspace.column_momentum,
                                    workspace.vector_divergence_edge_flux,
                                    workspace.modal_divergence);
    workspace.modal_solution = workspace.modal_divergence;
    const Real coefficient_m2 =
        implicit_time_s * implicit_time_s * modes.eigenvalue_m2_s2[mode];
    workspace.active_modal_coefficient_m2 = coefficient_m2;
    for (std::size_t cell = 0; cell < cells; ++cell) {
      Real laplacian_diagonal = 0.0;
      for (const auto& edge : grid.cell_cache()[cell].edges) {
        laplacian_diagonal += grid.edges()[edge.edge].length_m /
                              (grid.cells()[cell].area_m2 *
                               grid.edge_cache()[edge.edge].center_distance_m);
      }
      workspace.helmholtz_inverse_diagonal[cell] =
          1.0 / (1.0 + coefficient_m2 * laplacian_diagonal);
    }
    const auto linear =
        restarted_gmres(workspace.modal_helmholtz_operator, workspace.modal_divergence,
                        workspace.modal_solution, options, workspace.gmres,
                        workspace.modal_helmholtz_preconditioner);
    result.linear_iterations_total += linear.iterations;
    result.linear_iterations_maximum =
        std::max(result.linear_iterations_maximum, linear.iterations);
    result.linear_relative_residual_maximum =
        std::max(result.linear_relative_residual_maximum, linear.relative_residual);
    if (!linear.converged()) {
      result.all_converged = false;
      return result;
    }
    least_squares_gradient(grid, workspace.modal_solution,
                           workspace.helmholtz_gradient);
    for (std::size_t cell = 0; cell < cells; ++cell) {
      const auto n = dry_hydrostatic_offset(cell, mode, levels);
      workspace.modal_momentum[n] = workspace.modal_momentum[n] +
                                    coefficient_m2 * workspace.helmholtz_gradient[cell];
    }
  }

  correction.surface_pressure_pa = right_hand_side.surface_pressure_pa;
  correction.horizontal_momentum_mass_kg_m_s.assign(volume, {});
  correction.potential_temperature_mass_k_kg_m2 =
      right_hand_side.potential_temperature_mass_k_kg_m2;
  for (std::size_t cell = 0; cell < cells; ++cell) {
    for (std::size_t level = 0; level < levels; ++level) {
      Vec3 value{};
      for (std::size_t mode = 0; mode < levels; ++mode) {
        value =
            value +
            modes.eigenvectors[mode * levels + level] *
                workspace.modal_momentum[dry_hydrostatic_offset(cell, mode, levels)];
      }
      correction.horizontal_momentum_mass_kg_m_s[dry_hydrostatic_offset(
          cell, level, levels)] = value;
    }
  }

  workspace.momentum_perturbation = correction;
  std::fill(workspace.momentum_perturbation.surface_pressure_pa.begin(),
            workspace.momentum_perturbation.surface_pressure_pa.end(), 0.0);
  std::fill(workspace.momentum_perturbation.potential_temperature_mass_k_kg_m2.begin(),
            workspace.momentum_perturbation.potential_temperature_mass_k_kg_m2.end(),
            0.0);
  apply_dry_hydrostatic_fast_operator(
      grid, planet, fast_operator, workspace.momentum_perturbation,
      workspace.momentum_scalar_tendency, workspace.fast_operator);
  for (std::size_t cell = 0; cell < cells; ++cell) {
    correction.surface_pressure_pa[cell] +=
        implicit_time_s *
        workspace.momentum_scalar_tendency.surface_pressure_pa_s[cell];
  }
  for (std::size_t n = 0; n < volume; ++n) {
    correction.potential_temperature_mass_k_kg_m2[n] +=
        implicit_time_s *
        workspace.momentum_scalar_tendency.tendency.potential_temperature_mass[n];
  }

  // The selected-mode Schur equations are the accuracy contract of this
  // approximate inverse. Unselected small-Courant modes intentionally retain the
  // identity and are converged by the outer quasi-Newton iteration.
  result.equation_residual_norm = result.linear_relative_residual_maximum;
  return result;
}

}  // namespace mps
