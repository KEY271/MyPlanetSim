#include "myplanetsim/dynamics/dry_hydrostatic_driver.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

#include "myplanetsim/dynamics/dry_hydrostatic_benchmarks.hpp"
#include "myplanetsim/dynamics/dry_hydrostatic_diffusion.hpp"
#include "myplanetsim/dynamics/dry_hydrostatic_flux.hpp"
#include "myplanetsim/dynamics/dry_hydrostatic_reconstruction.hpp"
#include "myplanetsim/dynamics/dry_hydrostatic_sources.hpp"
#include "myplanetsim/numerics/additive_steppers.hpp"
#include "myplanetsim/numerics/diffusion_stability.hpp"
#include "myplanetsim/physics/dry_convective_adjustment.hpp"
#include "myplanetsim/physics/gray_radiation_coupling.hpp"
#include "myplanetsim/physics/held_suarez.hpp"
#include "myplanetsim/physics/moist_thermodynamics.hpp"
#include "myplanetsim/physics/physics_schedule.hpp"
#include "myplanetsim/physics/planetary_newtonian.hpp"
#include "myplanetsim/physics/surface_energy_balance.hpp"
namespace mps {
namespace {

[[nodiscard]] DryHydrostaticState euler_update(const DryHydrostaticState& base,
                                               const DryHydrostaticRhs& tendency,
                                               const Real scale) {
  DryHydrostaticState result = base;
  for (std::size_t cell = 0; cell < result.surface_pressure_pa.size(); ++cell) {
    result.surface_pressure_pa[cell] += scale * tendency.surface_pressure_pa_s[cell];
  }
  for (std::size_t n = 0; n < result.horizontal_momentum_mass_kg_m_s.size(); ++n) {
    result.horizontal_momentum_mass_kg_m_s[n] =
        result.horizontal_momentum_mass_kg_m_s[n] +
        scale * tendency.tendency.momentum[n];
    result.potential_temperature_mass_k_kg_m2[n] +=
        scale * tendency.tendency.potential_temperature_mass[n];
  }
  for (std::size_t n = 0; n < result.tracer_mass_kg_m2.size(); ++n)
    result.tracer_mass_kg_m2[n] += scale * tendency.tendency.tracer_mass[n];
  for (std::size_t cell = 0; cell < result.surface_temperature_k.size(); ++cell)
    result.surface_temperature_k[cell] +=
        scale * tendency.surface_temperature_k_s[cell];
  return result;
}

void add_scaled_rhs(DryHydrostaticState& state, const DryHydrostaticRhs& tendency,
                    const Real scale) {
  for (std::size_t cell = 0; cell < state.surface_pressure_pa.size(); ++cell)
    state.surface_pressure_pa[cell] += scale * tendency.surface_pressure_pa_s[cell];
  for (std::size_t n = 0; n < state.horizontal_momentum_mass_kg_m_s.size(); ++n) {
    state.horizontal_momentum_mass_kg_m_s[n] =
        state.horizontal_momentum_mass_kg_m_s[n] +
        scale * tendency.tendency.momentum[n];
    state.potential_temperature_mass_k_kg_m2[n] +=
        scale * tendency.tendency.potential_temperature_mass[n];
  }
  for (std::size_t n = 0; n < state.tracer_mass_kg_m2.size(); ++n)
    state.tracer_mass_kg_m2[n] += scale * tendency.tendency.tracer_mass[n];
  // The surface temperature is prognostic only when a surface model is active.
  for (std::size_t cell = 0; cell < state.surface_temperature_k.size(); ++cell)
    state.surface_temperature_k[cell] += scale * tendency.surface_temperature_k_s[cell];
}

void add_scaled_fast_tendency(DryHydrostaticState& state,
                              const DryHydrostaticFastTendency& tendency,
                              const Real scale) {
  for (std::size_t cell = 0; cell < state.surface_pressure_pa.size(); ++cell)
    state.surface_pressure_pa[cell] += scale * tendency.surface_pressure_pa_s[cell];
  for (std::size_t n = 0; n < state.horizontal_momentum_mass_kg_m_s.size(); ++n) {
    state.horizontal_momentum_mass_kg_m_s[n] =
        state.horizontal_momentum_mass_kg_m_s[n] +
        scale * tendency.tendency.momentum[n];
    state.potential_temperature_mass_k_kg_m2[n] +=
        scale * tendency.tendency.potential_temperature_mass[n];
  }
}

[[nodiscard]] DryHydrostaticState ark2_stage_right_hand_side(
    const DryHydrostaticState& initial,
    const std::array<const DryHydrostaticRhs*, 3>& full_tendencies,
    const std::array<const DryHydrostaticFastTendency*, 3>& fast_tendencies,
    const std::size_t stage, const Real time_step_s) {
  if (stage == 0 || stage >= 3)
    throw std::invalid_argument("ARK2 stage index must be one or two");
  const auto table = ark2_imex_table();
  DryHydrostaticState result = initial;
  for (std::size_t previous = 0; previous < stage; ++previous) {
    if (full_tendencies[previous] == nullptr || fast_tendencies[previous] == nullptr)
      throw std::invalid_argument("ARK2 stage tendency is unavailable");
    add_scaled_rhs(result, *full_tendencies[previous],
                   time_step_s * table.explicit_matrix[stage][previous]);
    add_scaled_fast_tendency(result, *fast_tendencies[previous],
                             time_step_s * (table.implicit_matrix[stage][previous] -
                                            table.explicit_matrix[stage][previous]));
  }
  return result;
}

[[nodiscard]] DryHydrostaticState ark2_final_update(
    const DryHydrostaticState& initial,
    const std::array<const DryHydrostaticRhs*, 3>& full_tendencies,
    const Real time_step_s) {
  const auto table = ark2_imex_table();
  DryHydrostaticState result = initial;
  for (std::size_t stage = 0; stage < 3; ++stage) {
    if (full_tendencies[stage] == nullptr)
      throw std::invalid_argument("ARK2 final tendency is unavailable");
    // The ARK2 explicit and implicit weights are identical, so N + L
    // recombines to the full nonlinear tendency at the final update.
    add_scaled_rhs(result, *full_tendencies[stage],
                   time_step_s * table.explicit_weights[stage]);
  }
  return result;
}

void assign_fast_solution(DryHydrostaticState& state,
                          const DryHydrostaticReferenceColumn& reference,
                          const DryHydrostaticFastPerturbation& solution) {
  const auto levels = reference.geometry.air_mass_kg_m2.size();
  for (std::size_t cell = 0; cell < state.surface_pressure_pa.size(); ++cell) {
    state.surface_pressure_pa[cell] =
        reference.surface_pressure_pa + solution.surface_pressure_pa[cell];
    for (std::size_t level = 0; level < levels; ++level) {
      const auto n = dry_hydrostatic_offset(cell, level, levels);
      state.horizontal_momentum_mass_kg_m_s[n] =
          solution.horizontal_momentum_mass_kg_m_s[n];
      state.potential_temperature_mass_k_kg_m2[n] =
          reference.potential_temperature_mass_k_kg_m2[level] +
          solution.potential_temperature_mass_k_kg_m2[n];
    }
  }
}

[[nodiscard]] DryHydrostaticState convex_update(const DryHydrostaticState& initial,
                                                const Real initial_weight,
                                                const DryHydrostaticState& stage,
                                                const DryHydrostaticRhs& tendency,
                                                const Real stage_weight,
                                                const Real time_step_s) {
  DryHydrostaticState result = initial;
  for (std::size_t cell = 0; cell < result.surface_pressure_pa.size(); ++cell) {
    result.surface_pressure_pa[cell] =
        initial_weight * initial.surface_pressure_pa[cell] +
        stage_weight * (stage.surface_pressure_pa[cell] +
                        time_step_s * tendency.surface_pressure_pa_s[cell]);
  }
  for (std::size_t n = 0; n < result.horizontal_momentum_mass_kg_m_s.size(); ++n) {
    result.horizontal_momentum_mass_kg_m_s[n] =
        initial_weight * initial.horizontal_momentum_mass_kg_m_s[n] +
        stage_weight * (stage.horizontal_momentum_mass_kg_m_s[n] +
                        time_step_s * tendency.tendency.momentum[n]);
    result.potential_temperature_mass_k_kg_m2[n] =
        initial_weight * initial.potential_temperature_mass_k_kg_m2[n] +
        stage_weight * (stage.potential_temperature_mass_k_kg_m2[n] +
                        time_step_s * tendency.tendency.potential_temperature_mass[n]);
  }
  for (std::size_t n = 0; n < result.tracer_mass_kg_m2.size(); ++n)
    result.tracer_mass_kg_m2[n] =
        initial_weight * initial.tracer_mass_kg_m2[n] +
        stage_weight * (stage.tracer_mass_kg_m2[n] +
                        time_step_s * tendency.tendency.tracer_mass[n]);
  for (std::size_t cell = 0; cell < result.surface_temperature_k.size(); ++cell)
    result.surface_temperature_k[cell] =
        initial_weight * initial.surface_temperature_k[cell] +
        stage_weight * (stage.surface_temperature_k[cell] +
                        time_step_s * tendency.surface_temperature_k_s[cell]);
  return result;
}

[[nodiscard]] bool pressure_is_in_range(const ExperimentConfig& config,
                                        const DryHydrostaticState& state) {
  return std::ranges::all_of(state.surface_pressure_pa, [&](const Real pressure) {
    return std::isfinite(pressure) &&
           pressure >= config.vertical.minimum_surface_pressure_pa &&
           pressure <= config.vertical.maximum_surface_pressure_pa;
  });
}

[[nodiscard]] Real pressure_limited_time_step(const ExperimentConfig& config,
                                              const DryHydrostaticState& state,
                                              const DryHydrostaticRhs& tendency,
                                              const Real maximum_time_step_s) {
  Real result = maximum_time_step_s;
  for (std::size_t cell = 0; cell < state.surface_pressure_pa.size(); ++cell) {
    const Real rate = tendency.surface_pressure_pa_s[cell];
    if (rate == 0.0) {
      continue;
    }
    const Real available = rate > 0.0 ? config.vertical.maximum_surface_pressure_pa -
                                            state.surface_pressure_pa[cell]
                                      : state.surface_pressure_pa[cell] -
                                            config.vertical.minimum_surface_pressure_pa;
    if (!(available > 0.0)) {
      throw std::runtime_error(
          "dry dynamics drives surface pressure outside the configured range");
    }
    const Real limit = available / std::abs(rate);
    if (limit < result) {
      result = std::nextafter(limit, 0.0);
    }
  }
  return result;
}

void halve_time_step(const DryHydrostaticState& state, Real& time_step_s) {
  time_step_s *= 0.5;
  if (!(time_step_s > 0.0) || state.time_s + time_step_s == state.time_s) {
    throw std::runtime_error(
        "dry stage pressure or CFL constraint produced a zero time step");
  }
}

[[nodiscard]] Real stable_time_step(const DryHydrostaticRhs& rhs) {
  return std::min({rhs.horizontal_stable_time_step_s, rhs.vertical_stable_time_step_s,
                   rhs.surface_stable_time_step_s, rhs.radiation_stable_time_step_s});
}

[[nodiscard]] DryConvectionDiagnostics apply_dry_convection(
    const ExperimentConfig& config, const CubedSphereGrid& grid,
    DryHydrostaticState& state, const DryHydrostaticDerived& derived,
    DryConvectiveAdjustmentResult& column_result,
    DryConvectiveAdjustmentWorkspace& column_workspace) {
  if (config.convection.kind == ConvectionKind::kNone) return {};
  const std::size_t cells = derived.cells;
  const std::size_t levels = derived.levels;
  if (cells != grid.cell_count() || levels == 0)
    throw std::invalid_argument("dry convection coupling shape mismatch");
  DryConvectionDiagnostics result{
      .minimum_theta_difference_before_k =
          levels > 1 ? std::numeric_limits<Real>::infinity() : 0.0,
      .minimum_theta_difference_after_k =
          levels > 1 ? std::numeric_limits<Real>::infinity() : 0.0};
  Real total_area = 0.0;
  for (std::size_t cell = 0; cell < cells; ++cell) {
    const std::size_t begin = cell * levels;
    const std::size_t half_begin = cell * (levels + 1);
    const DryConvectiveAdjustmentInput input{
        .potential_temperature_k =
            std::span<const Real>(derived.potential_temperature_k)
                .subspan(begin, levels),
        .air_mass_kg_m2 =
            std::span<const Real>(derived.air_mass_kg_m2).subspan(begin, levels),
        .exner_full = std::span<const Real>(derived.exner_full).subspan(begin, levels),
        .exner_half =
            std::span<const Real>(derived.exner_half).subspan(half_begin, levels + 1),
        .heat_capacity_cp_j_kg_k = config.planet.heat_capacity_cp_j_kg_k,
        .heat_capacity_cv_j_kg_k = config.planet.heat_capacity_cv_j_kg_k(),
        .stability_tolerance_k = config.convection.stability_tolerance_k};
    dry_convective_adjustment(input, column_result, column_workspace);
    const auto& column = column_result.diagnostics;
    for (std::size_t level = 0; level < levels; ++level)
      state.potential_temperature_mass_k_kg_m2[begin + level] =
          derived.air_mass_kg_m2[begin + level] *
          column_result.adjusted_potential_temperature_k[level];
    const Real area = grid.cells()[cell].area_m2;
    total_area += area;
    result.minimum_theta_difference_before_k =
        std::min(result.minimum_theta_difference_before_k,
                 column.minimum_theta_difference_before_k);
    result.minimum_theta_difference_after_k =
        std::min(result.minimum_theta_difference_after_k,
                 column.minimum_theta_difference_after_k);
    result.unstable_interface_fraction_before +=
        area * column.unstable_interface_fraction_before;
    result.unstable_interface_fraction_after +=
        area * column.unstable_interface_fraction_after;
    if (column.adjusted_layer_count > 0) ++result.adjusted_column_count;
    result.adjusted_layer_count += column.adjusted_layer_count;
    result.adjusted_block_count += column.adjusted_block_count;
    result.maximum_temperature_increment_k = std::max(
        result.maximum_temperature_increment_k, column.maximum_temperature_increment_k);
    result.enthalpy_change_j += area * column.enthalpy_change_j_m2;
    result.dry_energy_attributed_change_j +=
        area * column.dry_energy_attributed_change_j_m2;
  }
  if (!(total_area > 0.0))
    throw std::runtime_error("dry convection grid area is invalid");
  result.unstable_interface_fraction_before /= total_area;
  result.unstable_interface_fraction_after /= total_area;
  return result;
}

void accumulate_boundary_layer_diagnostics(DryMixingStepDiagnostics& total,
                                           const DryMixingStepDiagnostics& local,
                                           const Real time_weight) {
  total.sensible_to_atmosphere_energy_j += local.sensible_to_atmosphere_energy_j;
  total.surface_stress_impulse_magnitude_n_s +=
      local.surface_stress_impulse_magnitude_n_s;
  total.mean_boundary_layer_height_m +=
      time_weight * local.mean_boundary_layer_height_m;
  total.maximum_momentum_diffusivity_m2_s = std::max(
      total.maximum_momentum_diffusivity_m2_s, local.maximum_momentum_diffusivity_m2_s);
  total.maximum_heat_diffusivity_m2_s = std::max(total.maximum_heat_diffusivity_m2_s,
                                                 local.maximum_heat_diffusivity_m2_s);
  total.shallow_unresolved_area_fraction +=
      time_weight * local.shallow_unresolved_area_fraction;
  total.model_top_area_fraction += time_weight * local.model_top_area_fraction;
  total.atmospheric_heat_change_j += local.atmospheric_heat_change_j;
  total.surface_heat_change_j += local.surface_heat_change_j;
  total.physical_shear_dissipation_j += local.physical_shear_dissipation_j;
  total.physical_surface_drag_dissipation_j +=
      local.physical_surface_drag_dissipation_j;
  total.backward_euler_dissipation_j += local.backward_euler_dissipation_j;
  total.returned_dissipation_heat_j += local.returned_dissipation_heat_j;
  total.heat_budget_residual_j += local.heat_budget_residual_j;
  total.kinetic_energy_change_j += local.kinetic_energy_change_j;
  total.kinetic_energy_identity_residual_j += local.kinetic_energy_identity_residual_j;
  total.momentum_budget_residual_n_s += local.momentum_budget_residual_n_s;
  total.tracer_mass_change_kg += local.tracer_mass_change_kg;
  total.evaporation_kg += local.evaporation_kg;
  total.runoff_kg += local.runoff_kg;
  total.ocean_water_change_kg += local.ocean_water_change_kg;
  total.external_outflow_kg += local.external_outflow_kg;
  total.moist_enthalpy_budget_residual_j += local.moist_enthalpy_budget_residual_j;
  total.water_budget_residual_kg += local.water_budget_residual_kg;
}

void accumulate_convection_diagnostics(DryConvectionDiagnostics& total,
                                       const DryConvectionDiagnostics& local,
                                       const Real time_weight, bool& initialized) {
  if (!initialized) {
    total.minimum_theta_difference_before_k = local.minimum_theta_difference_before_k;
    total.minimum_theta_difference_after_k = local.minimum_theta_difference_after_k;
    initialized = true;
  } else {
    total.minimum_theta_difference_before_k =
        std::min(total.minimum_theta_difference_before_k,
                 local.minimum_theta_difference_before_k);
    total.minimum_theta_difference_after_k = std::min(
        total.minimum_theta_difference_after_k, local.minimum_theta_difference_after_k);
  }
  total.unstable_interface_fraction_before +=
      time_weight * local.unstable_interface_fraction_before;
  total.unstable_interface_fraction_after +=
      time_weight * local.unstable_interface_fraction_after;
  total.adjusted_column_count += local.adjusted_column_count;
  total.adjusted_layer_count += local.adjusted_layer_count;
  total.adjusted_block_count += local.adjusted_block_count;
  total.maximum_temperature_increment_k = std::max(
      total.maximum_temperature_increment_k, local.maximum_temperature_increment_k);
  total.enthalpy_change_j += local.enthalpy_change_j;
  total.dry_energy_attributed_change_j += local.dry_energy_attributed_change_j;
}

void accumulate_moist_physics_diagnostics(MoistPhysicsStepDiagnostics& total,
                                          const MoistPhysicsStepDiagnostics& local) {
  total.evaporation_kg += local.evaporation_kg;
  total.convective_precipitation_kg += local.convective_precipitation_kg;
  total.grid_scale_precipitation_kg += local.grid_scale_precipitation_kg;
  total.runoff_kg += local.runoff_kg;
  total.ocean_water_change_kg += local.ocean_water_change_kg;
  total.external_outflow_kg += local.external_outflow_kg;
  total.water_budget_residual_kg += local.water_budget_residual_kg;
  total.moist_enthalpy_budget_residual_j += local.moist_enthalpy_budget_residual_j;
  total.maximum_relative_humidity =
      std::max(total.maximum_relative_humidity, local.maximum_relative_humidity);
  total.maximum_temperature_increment_k = std::max(
      total.maximum_temperature_increment_k, local.maximum_temperature_increment_k);
  total.maximum_vapor_increment =
      std::max(total.maximum_vapor_increment, local.maximum_vapor_increment);
  total.cape_area_time_integral_j_m2_s_kg += local.cape_area_time_integral_j_m2_s_kg;
  total.cin_area_time_integral_j_m2_s_kg += local.cin_area_time_integral_j_m2_s_kg;
  total.lcl_pressure_area_time_integral_pa_m2_s +=
      local.lcl_pressure_area_time_integral_pa_m2_s;
  total.convection_top_pressure_area_time_integral_pa_m2_s +=
      local.convection_top_pressure_area_time_integral_pa_m2_s;
  total.active_area_time_m2_s += local.active_area_time_m2_s;
  total.deep_area_time_m2_s += local.deep_area_time_m2_s;
  total.shallow_area_time_m2_s += local.shallow_area_time_m2_s;
  total.inactive_area_time_m2_s += local.inactive_area_time_m2_s;
  total.model_top_area_time_m2_s += local.model_top_area_time_m2_s;
  total.deep_column_count += local.deep_column_count;
  total.shallow_column_count += local.shallow_column_count;
  total.inactive_column_count += local.inactive_column_count;
  total.sbm_diagnostic_column_count += local.sbm_diagnostic_column_count;
  total.sbm_relaxation_column_count += local.sbm_relaxation_column_count;
}

void resize_zero_rhs_term(DryHydrostaticRhsTerm& term, const std::size_t cells,
                          const std::size_t levels, const std::size_t tracer_count) {
  const auto cell_levels = cells * levels;
  term.surface_pressure_pa_s.assign(cells, 0.0);
  term.tendency.air_mass.assign(cell_levels, 0.0);
  term.tendency.momentum.assign(cell_levels, {});
  term.tendency.potential_temperature_mass.assign(cell_levels, 0.0);
  term.tendency.tracer_mass.assign(tracer_count * cell_levels, 0.0);
  term.surface_temperature_k_s.assign(cells, 0.0);
}

void resize_zero_rhs_components(DryHydrostaticRhsComponents& components,
                                const std::size_t cells, const std::size_t levels,
                                const std::size_t tracer_count) {
  resize_zero_rhs_term(components.horizontal_transport, cells, levels, tracer_count);
  resize_zero_rhs_term(components.vertical_transport, cells, levels, tracer_count);
  resize_zero_rhs_term(components.pressure_gradient, cells, levels, tracer_count);
  resize_zero_rhs_term(components.coriolis, cells, levels, tracer_count);
  resize_zero_rhs_term(components.diffusion, cells, levels, tracer_count);
  resize_zero_rhs_term(components.physics, cells, levels, tracer_count);
}

[[nodiscard]] bool surface_temperature_is_positive(const DryHydrostaticState& state) {
  return std::ranges::all_of(state.surface_temperature_k, [](const Real value) {
    return value > 0.0 && std::isfinite(value);
  });
}

[[nodiscard]] DryHydrostaticState crank_nicolson_residual(
    const DryHydrostaticState& initial, const DryHydrostaticState& candidate,
    const DryHydrostaticRhs& initial_rhs, const DryHydrostaticRhs& candidate_rhs,
    const Real time_step_s, const Real implicit_weight,
    const std::span<const Vec3> cell_centres, const std::size_t levels) {
  DryHydrostaticState residual = candidate;
  const Real initial_weight = 1.0 - implicit_weight;
  for (std::size_t cell = 0; cell < residual.surface_pressure_pa.size(); ++cell) {
    residual.surface_pressure_pa[cell] =
        candidate.surface_pressure_pa[cell] - initial.surface_pressure_pa[cell] -
        time_step_s * (initial_weight * initial_rhs.surface_pressure_pa_s[cell] +
                       implicit_weight * candidate_rhs.surface_pressure_pa_s[cell]);
  }
  for (std::size_t n = 0; n < residual.horizontal_momentum_mass_kg_m_s.size(); ++n) {
    residual.horizontal_momentum_mass_kg_m_s[n] =
        candidate.horizontal_momentum_mass_kg_m_s[n] -
        initial.horizontal_momentum_mass_kg_m_s[n] -
        time_step_s * (initial_weight * initial_rhs.tendency.momentum[n] +
                       implicit_weight * candidate_rhs.tendency.momentum[n]);
    residual.horizontal_momentum_mass_kg_m_s[n] = project_tangent(
        residual.horizontal_momentum_mass_kg_m_s[n], cell_centres[n / levels]);
    residual.potential_temperature_mass_k_kg_m2[n] =
        candidate.potential_temperature_mass_k_kg_m2[n] -
        initial.potential_temperature_mass_k_kg_m2[n] -
        time_step_s *
            (initial_weight * initial_rhs.tendency.potential_temperature_mass[n] +
             implicit_weight * candidate_rhs.tendency.potential_temperature_mass[n]);
  }
  for (std::size_t n = 0; n < residual.tracer_mass_kg_m2.size(); ++n)
    residual.tracer_mass_kg_m2[n] =
        candidate.tracer_mass_kg_m2[n] - initial.tracer_mass_kg_m2[n] -
        time_step_s * (initial_weight * initial_rhs.tendency.tracer_mass[n] +
                       implicit_weight * candidate_rhs.tendency.tracer_mass[n]);
  for (std::size_t cell = 0; cell < residual.surface_temperature_k.size(); ++cell) {
    residual.surface_temperature_k[cell] =
        candidate.surface_temperature_k[cell] - initial.surface_temperature_k[cell] -
        time_step_s * (initial_weight * initial_rhs.surface_temperature_k_s[cell] +
                       implicit_weight * candidate_rhs.surface_temperature_k_s[cell]);
  }
  return residual;
}

struct ScaledCrankNicolsonResidual {
  Real norm = 0.0;
  Real rms_norm = 0.0;
  const char* component = "none";
};

[[nodiscard]] ScaledCrankNicolsonResidual scaled_crank_nicolson_residual(
    const DryHydrostaticState& residual, const DryHydrostaticState& initial,
    const DryHydrostaticReferenceColumn& reference,
    const DryHydrostaticExternalModeOperator& external) {
  ScaledCrankNicolsonResidual result;
  Real sum_of_squares = 0.0;
  std::size_t sample_count = 0;
  const auto consider = [&](const Real value, const char* const component) {
    if (value > result.norm) {
      result.norm = value;
      result.component = component;
    }
    sum_of_squares += value * value;
    ++sample_count;
  };
  for (std::size_t cell = 0; cell < residual.surface_pressure_pa.size(); ++cell) {
    consider(std::abs(residual.surface_pressure_pa[cell]) /
                 std::max({1.0, reference.surface_pressure_pa,
                           std::abs(initial.surface_pressure_pa[cell])}),
             "surface_pressure");
  }
  const auto levels = reference.geometry.air_mass_kg_m2.size();
  for (std::size_t n = 0; n < residual.horizontal_momentum_mass_kg_m_s.size(); ++n) {
    const auto level = n % levels;
    const Real momentum_scale = std::max(
        1.0, reference.geometry.air_mass_kg_m2[level] * external.phase_speed_m_s);
    consider(
        norm(residual.horizontal_momentum_mass_kg_m_s[n]) /
            std::max(momentum_scale, norm(initial.horizontal_momentum_mass_kg_m_s[n])),
        "momentum");
    consider(std::abs(residual.potential_temperature_mass_k_kg_m2[n]) /
                 std::max({1.0, reference.potential_temperature_mass_k_kg_m2[level],
                           std::abs(initial.potential_temperature_mass_k_kg_m2[n])}),
             "potential_temperature_mass");
  }
  for (std::size_t n = 0; n < residual.tracer_mass_kg_m2.size(); ++n)
    consider(std::abs(residual.tracer_mass_kg_m2[n]) /
                 std::max(1.0, std::abs(initial.tracer_mass_kg_m2[n])),
             "tracer_mass");
  for (std::size_t cell = 0; cell < residual.surface_temperature_k.size(); ++cell) {
    consider(std::abs(residual.surface_temperature_k[cell]) /
                 std::max({1.0, reference.temperature_k,
                           std::abs(initial.surface_temperature_k[cell])}),
             "surface_temperature");
  }
  result.rms_norm = std::sqrt(sum_of_squares / static_cast<Real>(sample_count));
  return result;
}

[[nodiscard]] DryHydrostaticFastPerturbation negative_fast_residual(
    const DryHydrostaticState& residual) {
  DryHydrostaticFastPerturbation result{
      .surface_pressure_pa = residual.surface_pressure_pa,
      .horizontal_momentum_mass_kg_m_s = residual.horizontal_momentum_mass_kg_m_s,
      .potential_temperature_mass_k_kg_m2 =
          residual.potential_temperature_mass_k_kg_m2};
  for (auto& value : result.surface_pressure_pa) value = -value;
  for (auto& value : result.horizontal_momentum_mass_kg_m_s) value = -value;
  for (auto& value : result.potential_temperature_mass_k_kg_m2) value = -value;
  return result;
}

void add_semi_implicit_correction(DryHydrostaticState& state,
                                  const DryHydrostaticFastPerturbation& fast_correction,
                                  const DryHydrostaticState& residual) {
  for (std::size_t cell = 0; cell < state.surface_pressure_pa.size(); ++cell)
    state.surface_pressure_pa[cell] += fast_correction.surface_pressure_pa[cell];
  for (std::size_t n = 0; n < state.horizontal_momentum_mass_kg_m_s.size(); ++n) {
    state.horizontal_momentum_mass_kg_m_s[n] =
        state.horizontal_momentum_mass_kg_m_s[n] +
        fast_correction.horizontal_momentum_mass_kg_m_s[n];
    state.potential_temperature_mass_k_kg_m2[n] +=
        fast_correction.potential_temperature_mass_k_kg_m2[n];
    // Tracer and surface reservoirs are outside L_ref, so their quasi-Newton
    // correction is the identity solve.
  }
  for (std::size_t n = 0; n < state.tracer_mass_kg_m2.size(); ++n)
    state.tracer_mass_kg_m2[n] -= residual.tracer_mass_kg_m2[n];
  for (std::size_t cell = 0; cell < state.surface_temperature_k.size(); ++cell)
    state.surface_temperature_k[cell] -= residual.surface_temperature_k[cell];
}

[[nodiscard]] Real semi_implicit_explicit_limit(const ExperimentConfig& config,
                                                const DryHydrostaticRhs& rhs) {
  const Real advective = rhs.horizontal_advective_stable_time_step_s *
                         config.dry_hydrostatic.advective_cfl /
                         config.dry_hydrostatic.cfl;
  return std::min({advective, rhs.diffusion_stable_time_step_s,
                   rhs.vertical_stable_time_step_s, rhs.surface_stable_time_step_s,
                   rhs.radiation_stable_time_step_s});
}

[[nodiscard]] Real courant_from_stable_time_step(const Real time_step_s,
                                                 const Real configured_cfl,
                                                 const Real stable_time_step_s) {
  if (std::isinf(stable_time_step_s)) return 0.0;
  if (!(stable_time_step_s > 0.0) || !std::isfinite(stable_time_step_s))
    throw std::runtime_error("stable time step is invalid for Courant diagnostics");
  return time_step_s * configured_cfl / stable_time_step_s;
}

[[nodiscard]] std::string nonlinear_failure_message(const Real initial_residual,
                                                    const Real residual,
                                                    const Real rms_residual,
                                                    const char* const component,
                                                    const Real time_step_s) {
  std::ostringstream message;
  message << std::setprecision(17)
          << "Crank-Nicolson iteration diverged: initial_residual=" << initial_residual
          << ", residual=" << residual << ", rms_residual=" << rms_residual
          << ", component=" << component << ", dt_s=" << time_step_s;
  return message.str();
}

[[nodiscard]] std::string modal_failure_message(const Real equation_residual) {
  std::ostringstream message;
  message << std::setprecision(17)
          << "modal linear solve did not converge: equation_residual="
          << equation_residual;
  return message.str();
}

[[nodiscard]] GrayRadiationDiagnostics weighted_radiation_diagnostics(
    const GrayRadiationDiagnostics& first, const Real first_weight,
    const GrayRadiationDiagnostics& second, const Real second_weight,
    const GrayRadiationDiagnostics& third = {}, const Real third_weight = 0.0) {
  const auto weighted = [&](const Real GrayRadiationDiagnostics::* member) {
    return first_weight * first.*member + second_weight * second.*member +
           third_weight * third.*member;
  };
  return {
      .toa_incoming_shortwave_power_w =
          weighted(&GrayRadiationDiagnostics::toa_incoming_shortwave_power_w),
      .toa_reflected_shortwave_power_w =
          weighted(&GrayRadiationDiagnostics::toa_reflected_shortwave_power_w),
      .toa_outgoing_longwave_power_w =
          weighted(&GrayRadiationDiagnostics::toa_outgoing_longwave_power_w),
      .toa_net_upward_power_w =
          weighted(&GrayRadiationDiagnostics::toa_net_upward_power_w),
      .surface_down_shortwave_power_w =
          weighted(&GrayRadiationDiagnostics::surface_down_shortwave_power_w),
      .surface_up_shortwave_power_w =
          weighted(&GrayRadiationDiagnostics::surface_up_shortwave_power_w),
      .surface_down_longwave_power_w =
          weighted(&GrayRadiationDiagnostics::surface_down_longwave_power_w),
      .surface_up_longwave_power_w =
          weighted(&GrayRadiationDiagnostics::surface_up_longwave_power_w),
      .atmospheric_shortwave_heating_power_w =
          weighted(&GrayRadiationDiagnostics::atmospheric_shortwave_heating_power_w),
      .atmospheric_longwave_heating_power_w =
          weighted(&GrayRadiationDiagnostics::atmospheric_longwave_heating_power_w),
      .surface_storage_rate_w =
          weighted(&GrayRadiationDiagnostics::surface_storage_rate_w),
      .sensible_to_atmosphere_power_w =
          weighted(&GrayRadiationDiagnostics::sensible_to_atmosphere_power_w),
      .internal_heat_power_w =
          weighted(&GrayRadiationDiagnostics::internal_heat_power_w),
      .interface_conservation_residual_w =
          weighted(&GrayRadiationDiagnostics::interface_conservation_residual_w),
      .dry_thermal_energy_rate_w =
          weighted(&GrayRadiationDiagnostics::dry_thermal_energy_rate_w),
      .rayleigh_drag_work_w = weighted(&GrayRadiationDiagnostics::rayleigh_drag_work_w),
  };
}

struct ScheduledRadiationStep {
  GrayRadiationDiagnostics rates{};
  Real stable_time_step_s = std::numeric_limits<Real>::infinity();
  std::vector<Real> surface_temperature_after_update;
  std::size_t column_call_count = 0;
  Real wall_seconds = 0.0;
};

}  // namespace

DryHydrostaticDriver::DryHydrostaticDriver(ExperimentConfig c,
                                           const std::size_t horizontal_tile_width)
    : config_(std::move(c)),
      grid_(config_.grid.cells_per_panel, config_.planet.radius_m),
      horizontal_tiles_(make_cubed_sphere_tiles(grid_, horizontal_tile_width,
                                                kDryHydrostaticTileHaloRadius)),
      coordinate_({config_.vertical.a_half_pa, config_.vertical.b_half},
                  config_.vertical.minimum_surface_pressure_pa,
                  config_.vertical.maximum_surface_pressure_pa,
                  config_.vertical.minimum_pressure_thickness_pa),
      orography_(make_surface_orography(config_.orography, grid_, config_.planet,
                                        config_.source_directory)) {
  if (config_.kind != ExperimentKind::kDryHydrostatic)
    throw std::invalid_argument("dry driver requires dry_hydrostatic config");
  if (config_.surface.has_value()) {
    surface_boundary_ =
        make_surface_boundary(*config_.surface, config_.orography, grid_,
                              config_.planet, config_.source_directory);
    orography_ = SurfaceOrography(
        std::vector<Real>(surface_boundary_->surface_geopotential_m2_s2().begin(),
                          surface_boundary_->surface_geopotential_m2_s2().end()),
        surface_boundary_->source_fingerprint());
  }
  // Reference-state well balancing is intentionally limited to the named hydrostatic
  // rest benchmark. Applying an initial-state correction to jets or forced cases would
  // silently modify their physical pressure force.
  if (config_.dry_hydrostatic.test_case == DryHydrostaticTestCase::kDcmip200Rest) {
    const auto reference_state = initial_state();
    pressure_reference_ = make_dry_hydrostatic_pressure_reference(
        grid_, diagnose(reference_state), config_.planet);
  }
  if (config_.dry_hydrostatic.time_integrator !=
      DryHydrostaticTimeIntegrator::kExplicitSspRk3) {
    if (!config_.semi_implicit.has_value())
      throw std::invalid_argument(
          "semi-implicit dry driver requires semi-implicit parameters");
    semi_implicit_reference_column_ = make_dry_hydrostatic_reference_column(
        coordinate_, config_.planet, *config_.semi_implicit);
    semi_implicit_fast_operator_ = make_dry_hydrostatic_fast_operator(
        coordinate_, config_.planet, *semi_implicit_reference_column_);
    semi_implicit_vertical_modes_ = make_dry_hydrostatic_vertical_modes(
        *semi_implicit_reference_column_, config_.planet,
        *semi_implicit_fast_operator_);
    semi_implicit_external_operator_ = make_dry_hydrostatic_external_mode_operator(
        *semi_implicit_reference_column_, *semi_implicit_vertical_modes_);
  }
}

void DryHydrostaticDriver::update_semi_implicit_reference(
    const DryHydrostaticState& state) const {
  // ADR 0016: keep the reference column horizontally uniform so a single vertical
  // structure eigendecomposition still serves every cell, but track the current state
  // so the reference-linear operator stays close to the true Jacobian.
  const auto levels = coordinate_.levels();
  const auto cells = grid_.cell_count();
  diagnose_dry_hydrostatic_state(state, coordinate_, config_.planet,
                                 orography_.surface_geopotential_m2_s2(),
                                 workspace_.derived, workspace_.diagnosis_workspace);
  std::vector<Real> profile(levels, 0.0);
  Real area = 0.0;
  Real surface_pressure = 0.0;
  for (std::size_t cell = 0; cell < cells; ++cell) {
    const Real cell_area = grid_.cells()[cell].area_m2;
    area += cell_area;
    surface_pressure += cell_area * state.surface_pressure_pa[cell];
    for (std::size_t level = 0; level < levels; ++level)
      profile[level] +=
          cell_area *
          workspace_.derived.temperature_k[dry_hydrostatic_offset(cell, level, levels)];
  }
  if (!(area > 0.0)) throw std::runtime_error("grid area is invalid");
  surface_pressure /= area;
  for (auto& value : profile) value /= area;
  semi_implicit_reference_column_ = make_dry_hydrostatic_reference_column(
      coordinate_, config_.planet, surface_pressure, profile);
  semi_implicit_fast_operator_ = make_dry_hydrostatic_fast_operator(
      coordinate_, config_.planet, *semi_implicit_reference_column_);
  semi_implicit_vertical_modes_ = make_dry_hydrostatic_vertical_modes(
      *semi_implicit_reference_column_, config_.planet, *semi_implicit_fast_operator_);
  semi_implicit_external_operator_ = make_dry_hydrostatic_external_mode_operator(
      *semi_implicit_reference_column_, *semi_implicit_vertical_modes_);
}
DryHydrostaticState DryHydrostaticDriver::initial_state() const {
  auto state =
      initialize_dry_hydrostatic_benchmark(config_, grid_, coordinate_, orography_);
  if (!config_.tracers.empty()) {
    const auto cells = grid_.cell_count();
    const auto levels = coordinate_.levels();
    const auto volume = cells * levels;
    const auto derived = diagnose(state);
    state.tracer_count = config_.tracers.size();
    state.tracer_mass_kg_m2.resize(state.tracer_count * volume);
    for (std::size_t tracer = 0; tracer < state.tracer_count; ++tracer)
      for (std::size_t n = 0; n < volume; ++n)
        state.tracer_mass_kg_m2[tracer * volume + n] =
            derived.air_mass_kg_m2[n] *
            (config_.tracers[tracer].initial_relative_humidity
                 ? *config_.tracers[tracer].initial_relative_humidity *
                       saturation_mixing_ratio(
                           derived.temperature_k[n], derived.pressure_pa[n],
                           {.gas_constant_dry_air_j_kg_k =
                                config_.planet.gas_constant_j_kg_k,
                            .heat_capacity_cp_j_kg_k =
                                config_.planet.heat_capacity_cp_j_kg_k})
                 : config_.tracers[tracer].initial_mixing_ratio);

    if (config_.moisture.kind == MoistureKind::kDiluteWater) {
      const TracerRegistry registry(config_.tracers);
      const auto water_vapor = *registry.water_vapor_index();
      const DiluteMoistThermodynamics moist_thermodynamics{
          .gas_constant_dry_air_j_kg_k = config_.planet.gas_constant_j_kg_k,
          .heat_capacity_cp_j_kg_k = config_.planet.heat_capacity_cp_j_kg_k};
      for (std::size_t n = 0; n < volume; ++n) {
        validate_dilute_moist_state(derived.temperature_k[n], derived.pressure_pa[n],
                                    state.tracer_mass_kg_m2[water_vapor * volume + n] /
                                        derived.air_mass_kg_m2[n],
                                    moist_thermodynamics);
      }
    }
  }
  if (config_.physics.kind == PhysicsKind::kSurfaceEnergyBalance ||
      config_.physics.kind == PhysicsKind::kGrayRadiation)
    state.surface_temperature_k.assign(grid_.cell_count(),
                                       config_.surface->initial_temperature_k);
  if (config_.moisture.kind == MoistureKind::kDiluteWater) {
    state.land_water_kg_m2.resize(grid_.cell_count());
    for (std::size_t cell = 0; cell < grid_.cell_count(); ++cell)
      if (surface_boundary_->land_fraction()[cell] > 0.0)
        state.land_water_kg_m2[cell] = config_.surface->hydrology_capacity_kg_m2 *
                                       config_.surface->hydrology_initial_fraction;
  }
  return state;
}
DryHydrostaticDerived DryHydrostaticDriver::diagnose(
    const DryHydrostaticState& s) const {
  return diagnose_dry_hydrostatic_state(s, coordinate_, config_.planet,
                                        orography_.surface_geopotential_m2_s2());
}
DryHydrostaticSources DryHydrostaticDriver::diagnose_sources(
    const DryHydrostaticDerived& d) const {
  DryHydrostaticSources result;
  DryHydrostaticSourcesWorkspace workspace;
  if (pressure_reference_.has_value())
    dry_hydrostatic_sources(grid_, d, config_.planet, *pressure_reference_, result,
                            workspace);
  else
    dry_hydrostatic_sources(grid_, d, config_.planet, result, workspace);
  return result;
}
DryHydrostaticRhs DryHydrostaticDriver::rhs(const DryHydrostaticState& s) const {
  DryHydrostaticRhs result;
  rhs(s, result);
  return result;
}

void DryHydrostaticDriver::rhs(const DryHydrostaticState& s,
                               DryHydrostaticRhs& result) const {
  rhs_with_components(s, result, nullptr);
}

DryHydrostaticRhsComponents DryHydrostaticDriver::rhs_components(
    const DryHydrostaticState& s) const {
  DryHydrostaticRhs result;
  DryHydrostaticRhsComponents components;
  rhs_with_components(s, result, &components);
  return components;
}

void DryHydrostaticDriver::rhs_with_components(
    const DryHydrostaticState& s, DryHydrostaticRhs& result,
    DryHydrostaticRhsComponents* const components,
    const bool compute_fast_wave_cfl) const {
  auto stamp = profile_enabled_ ? std::chrono::steady_clock::now()
                                : std::chrono::steady_clock::time_point{};
  const auto finish_region = [&](std::size_t region) {
    if (!profile_enabled_) return;
    const auto now = std::chrono::steady_clock::now();
    rhs_profile_.seconds[region] += std::chrono::duration<Real>(now - stamp).count();
    stamp = now;
  };
  if (profile_enabled_) {
    ++rhs_profile_.calls;
    std::size_t halo_references = 0;
    std::size_t index_references = 0;
    for (const auto& tile : horizontal_tiles_) {
      halo_references += tile.halo_cells.size();
      index_references +=
          tile.interior_cells.size() + tile.halo_cells.size() + tile.owned_edges.size();
    }
    rhs_profile_.tile_count = horizontal_tiles_.size();
    rhs_profile_.tile_halo_cell_references = halo_references;
    rhs_profile_.tile_metadata_bytes =
        horizontal_tiles_.size() * sizeof(CubedSphereTile) +
        index_references * sizeof(std::size_t);
  }
  auto& workspace = workspace_;
  diagnose_dry_hydrostatic_state(s, coordinate_, config_.planet,
                                 orography_.surface_geopotential_m2_s2(),
                                 workspace.derived, workspace.diagnosis_workspace);
  const auto& d = workspace.derived;
  auto C = d.cells, K = d.levels;
  if (components != nullptr)
    resize_zero_rhs_components(*components, C, K, d.tracer_count);
  auto& h = workspace.horizontal_tendency;
  h.air_mass.assign(C * K, 0.0);
  h.momentum.assign(C * K, {});
  h.potential_temperature_mass.assign(C * K, 0.0);
  h.tracer_mass.assign(d.tracer_count * C * K, 0.0);
  // Per-cell Courant condition (ADR 0011): dt * sum_f(lambda_f * L_f) / A_cell <= cfl,
  // the same definition the transport and shallow-water solvers already use. The
  // previous per-edge form was about four times weaker on a quadrilateral cell.
  if (compute_fast_wave_cfl) workspace.face_fast_wave_speed_length.assign(C * K, 0.0);
  workspace.face_advective_speed_length.assign(C * K, 0.0);
  auto& face_fast_wave_speed_length = workspace.face_fast_wave_speed_length;
  auto& face_advective_speed_length = workspace.face_advective_speed_length;
  finish_region(0);
  prepare_dry_hydrostatic_reconstruction(
      grid_, d, config_.dry_hydrostatic.reconstruction, config_.dry_hydrostatic.limiter,
      workspace.prepared_reconstruction, workspace.reconstruction_workspace,
      compute_fast_wave_cfl);
  if (profile_enabled_) {
    const auto scalar_bytes = [](const auto& field) {
      return field.gradient.size() * sizeof(Vec3) + field.factor.size() * sizeof(Real);
    };
    rhs_profile_.prepared_reconstruction_bytes = std::max(
        rhs_profile_.prepared_reconstruction_bytes,
        scalar_bytes(workspace.prepared_reconstruction.air_mass) +
            scalar_bytes(workspace.prepared_reconstruction.potential_temperature) +
            scalar_bytes(workspace.prepared_reconstruction.temperature) +
            scalar_bytes(workspace.prepared_reconstruction.tracer) +
            workspace.prepared_reconstruction.velocity_gradient.size() *
                sizeof(TangentVectorGradient) +
            workspace.prepared_reconstruction.velocity_factor.size() * sizeof(Real) +
            workspace.prepared_reconstruction.tracer_is_constant.size() *
                sizeof(std::uint8_t));
    rhs_profile_.eliminated_face_state_bytes =
        std::max(rhs_profile_.eliminated_face_state_bytes,
                 grid_.edge_count() * K * sizeof(DryHydrostaticFaceStates) +
                     2 * d.tracer_count * grid_.edge_count() * K * sizeof(Real));
  }
  finish_region(1);
  const auto& reconstructed = workspace.prepared_reconstruction;
  const auto edges = grid_.edges();
  const auto edge_count = edges.size();
  workspace.edge_flux.resize(edge_count * K);
  workspace.edge_tracer_flux.resize(d.tracer_count * edge_count * K);
  workspace.edge_flux_failures.resize(edge_count);
  if (profile_enabled_)
    rhs_profile_.edge_flux_bytes =
        std::max(rhs_profile_.edge_flux_bytes,
                 workspace.edge_flux.size() * sizeof(DryHydrostaticEdgeFlux) +
                     workspace.edge_tracer_flux.size() * sizeof(Real));
  std::fill(workspace.edge_flux_failures.begin(), workspace.edge_flux_failures.end(),
            nullptr);
#if defined(MPS_ENABLE_OPENMP)
#pragma omp parallel for schedule(static)
#endif
  for (std::ptrdiff_t tile_index = 0;
       tile_index < static_cast<std::ptrdiff_t>(horizontal_tiles_.size());
       ++tile_index) {
    for (const auto edge_id :
         horizontal_tiles_[static_cast<std::size_t>(tile_index)].owned_edges) {
      const auto& e = edges[edge_id];
      const auto& cached_edge = grid_.edge_cache()[e.id];
      const EdgeTangentBasis basis{cached_edge.normal, cached_edge.tangent};
      try {
        for (std::size_t k = 0; k < K; ++k) {
          const auto face =
              reconstruct_dry_hydrostatic_edge(grid_, d, reconstructed, e.id, k);
          workspace.edge_flux[e.id * K + k] = rusanov_dry_hydrostatic_flux(
              face.left, face.right, basis, config_.planet.gas_constant_j_kg_k,
              config_.planet.heat_capacity_cp_j_kg_k, compute_fast_wave_cfl);
          for (std::size_t tracer = 0; tracer < d.tracer_count; ++tracer) {
            const auto q = (tracer * edge_count + e.id) * K + k;
            workspace.edge_tracer_flux[q] = rusanov_dry_hydrostatic_tracer_flux(
                face.left, face.right,
                reconstruct_dry_hydrostatic_tracer_face(grid_, d, reconstructed, true,
                                                        tracer, e.id, k),
                reconstruct_dry_hydrostatic_tracer_face(grid_, d, reconstructed, false,
                                                        tracer, e.id, k),
                basis);
          }
        }
      } catch (...) {
        workspace.edge_flux_failures[edge_id] = std::current_exception();
      }
    }
  }
  for (std::size_t edge = 0; edge < edge_count; ++edge)
    if (workspace.edge_flux_failures[edge] != nullptr)
      std::rethrow_exception(workspace.edge_flux_failures[edge]);

#if defined(MPS_ENABLE_OPENMP)
#pragma omp parallel for schedule(static)
#endif
  for (std::ptrdiff_t tile_index = 0;
       tile_index < static_cast<std::ptrdiff_t>(horizontal_tiles_.size());
       ++tile_index) {
    for (const auto cell :
         horizontal_tiles_[static_cast<std::size_t>(tile_index)].interior_cells) {
      std::array<std::size_t, 4> cell_edges{};
      for (std::size_t side = 0; side < cell_edges.size(); ++side)
        cell_edges[side] = grid_.cell_cache()[cell].edges[side].edge;
      std::ranges::sort(cell_edges);
      for (const auto edge : cell_edges) {
        const auto& e = edges[edge];
        const auto& cached = grid_.edge_cache()[edge];
        const Real sign = cached.left_cell == cell ? -1.0 : 1.0;
        const Real scale = sign * e.length_m / grid_.cells()[cell].area_m2;
        for (std::size_t k = 0; k < K; ++k) {
          const auto n = dry_hydrostatic_offset(cell, k, K);
          const auto& flux = workspace.edge_flux[edge * K + k];
          h.air_mass[n] += scale * flux.air_mass_kg_m_s;
          h.momentum[n] = h.momentum[n] + scale * flux.momentum_kg_s2;
          h.potential_temperature_mass[n] +=
              scale * flux.potential_temperature_mass_k_kg_m_s;
          for (std::size_t tracer = 0; tracer < d.tracer_count; ++tracer) {
            const auto q = dry_hydrostatic_tracer_offset(tracer, cell, k, C, K);
            const auto edge_q = (tracer * edge_count + edge) * K + k;
            h.tracer_mass[q] += scale * workspace.edge_tracer_flux[edge_q];
          }
          if (compute_fast_wave_cfl)
            face_fast_wave_speed_length[n] += e.length_m * flux.maximum_wave_speed_m_s;
          face_advective_speed_length[n] +=
              e.length_m * flux.maximum_dissipation_speed_m_s;
        }
      }
    }
  }
  Real fast_wave_dt = std::numeric_limits<Real>::infinity();
  Real advective_dt = std::numeric_limits<Real>::infinity();
  for (std::size_t c = 0; c < C; ++c) {
    for (std::size_t k = 0; k < K; ++k) {
      const auto n = dry_hydrostatic_offset(c, k, K);
      if (compute_fast_wave_cfl) {
        const Real fast_wave_denominator = face_fast_wave_speed_length[n];
        if (fast_wave_denominator > 0.0)
          fast_wave_dt = std::min(fast_wave_dt, config_.dry_hydrostatic.cfl *
                                                    grid_.cells()[c].area_m2 /
                                                    fast_wave_denominator);
      }
      const Real advective_denominator = face_advective_speed_length[n];
      if (advective_denominator > 0.0)
        advective_dt = std::min(advective_dt, config_.dry_hydrostatic.cfl *
                                                  grid_.cells()[c].area_m2 /
                                                  advective_denominator);
    }
  }
  finish_region(2);
  couple_dry_hydrostatic_columns(
      s, d, h, coordinate_.coefficients().b_half, config_.planet.gravity_m_s2,
      config_.vertical.transport_scheme, config_.vertical.limiter,
      // Nonlinear minmod switches in theta transport degrade the cancellation between
      // thermodynamic transport and hydrostatic pressure work. Keep the configured
      // limiter for momentum and passive tracer, but use the smoother linear theta
      // reconstruction. Stage validation still rejects non-positive temperatures.
      VerticalLimiterKind::kNone, workspace.coupling, workspace.coupling_workspace);
  auto& coupled = workspace.coupling;
  if (components != nullptr) {
    components->horizontal_transport.tendency = h;
    components->vertical_transport.surface_pressure_pa_s =
        coupled.surface_pressure_pa_s;
    for (std::size_t n = 0; n < C * K; ++n) {
      components->vertical_transport.tendency.air_mass[n] =
          coupled.tendency.air_mass[n] - h.air_mass[n];
      components->vertical_transport.tendency.momentum[n] =
          coupled.tendency.momentum[n] - h.momentum[n];
      components->vertical_transport.tendency.potential_temperature_mass[n] =
          coupled.tendency.potential_temperature_mass[n] -
          h.potential_temperature_mass[n];
    }
    for (std::size_t n = 0; n < d.tracer_count * C * K; ++n)
      components->vertical_transport.tendency.tracer_mass[n] =
          coupled.tendency.tracer_mass[n] - h.tracer_mass[n];
  }
  Real vertical_dt = std::numeric_limits<Real>::infinity();
  for (std::size_t c = 0; c < C; ++c) {
    const auto begin = c * K;
    const auto interface_begin = c * (K + 1);
    const std::span<const Real> air_mass(d.air_mass_kg_m2.data() + begin, K);
    const std::span<const Real> interface_flux(
        coupled.interface_mass_flux_kg_m2_s.data() + interface_begin, K + 1);
    const std::span<const Real> horizontal_air_mass(h.air_mass.data() + begin, K);
    const Real unit_step_cfl =
        vertical_maximum_cfl(air_mass, interface_flux, horizontal_air_mass, 1.0);
    if (unit_step_cfl > 0.0)
      vertical_dt = std::min(vertical_dt, config_.vertical.cfl / unit_step_cfl);
  }
  finish_region(3);
  if (pressure_reference_.has_value())
    dry_hydrostatic_sources(grid_, d, config_.planet, *pressure_reference_,
                            workspace.sources, workspace.sources_workspace);
  else
    dry_hydrostatic_sources(grid_, d, config_.planet, workspace.sources,
                            workspace.sources_workspace);
  const auto& sources = workspace.sources;
  if (components != nullptr) {
    components->pressure_gradient.tendency.momentum = sources.pressure_gradient_kg_m_s2;
    components->coriolis.tendency.momentum = sources.coriolis_kg_m_s2;
  }
  for (std::size_t n = 0; n < C * K; ++n)
    coupled.tendency.momentum[n] = coupled.tendency.momentum[n] +
                                   sources.pressure_gradient_kg_m_s2[n] +
                                   sources.coriolis_kg_m_s2[n];
  finish_region(4);
  Real diffusion_rate = 0.0;
  Real diffusion_dt = std::numeric_limits<Real>::infinity();
  if (config_.dry_hydrostatic.diffusion_kind != DiffusionKind::kNone) {
    const auto diffusion = dry_hydrostatic_diffusion_tendency(
        grid_, d, config_.dry_hydrostatic.diffusion_kind,
        config_.dry_hydrostatic.diffusion_coefficient, config_.tracers);
    for (std::size_t n = 0; n < C * K; ++n) {
      coupled.tendency.momentum[n] =
          coupled.tendency.momentum[n] + diffusion.momentum[n];
      coupled.tendency.potential_temperature_mass[n] +=
          diffusion.potential_temperature_mass[n];
    }
    for (std::size_t n = 0; n < d.tracer_count * C * K; ++n)
      coupled.tendency.tracer_mass[n] += diffusion.tracer_mass[n];
    diffusion_rate = diffusion.kinetic_energy_rate_w;
    if (components != nullptr) {
      components->diffusion.tendency.momentum = diffusion.momentum;
      components->diffusion.tendency.potential_temperature_mass =
          diffusion.potential_temperature_mass;
      components->diffusion.tendency.tracer_mass = diffusion.tracer_mass;
    }
    diffusion_dt =
        stable_diffusion_time_step(grid_, config_.dry_hydrostatic.diffusion_kind,
                                   config_.dry_hydrostatic.diffusion_coefficient,
                                   std::numeric_limits<Real>::max());
  }
  finish_region(5);
  HeldSuarezDiagnostics physics_diagnostics{};
  SurfaceEnergyDiagnostics surface_diagnostics{};
  GrayRadiationDiagnostics radiation_diagnostics{};
  std::size_t radiation_column_call_count = 0;
  Real radiation_wall_seconds = 0.0;
  Real surface_dt = std::numeric_limits<Real>::infinity();
  Real radiation_dt = std::numeric_limits<Real>::infinity();
  // A split local-physics path can own a nonempty surface state while contributing
  // no surface tendency to the dynamics RHS.
  result.surface_temperature_k_s.assign(s.surface_temperature_k.size(), 0.0);
  if (config_.physics.kind == PhysicsKind::kHeldSuarez) {
    held_suarez_tendency(grid_, coordinate_, d, s.surface_pressure_pa, config_.planet,
                         workspace.atmospheric_physics,
                         workspace.atmospheric_physics_workspace);
    const auto& physics = workspace.atmospheric_physics;
    if (components != nullptr) {
      components->physics.tendency.momentum = physics.horizontal_momentum_mass_kg_m_s2;
      components->physics.tendency.potential_temperature_mass =
          physics.potential_temperature_mass_k_kg_m2_s;
    }
    for (std::size_t n = 0; n < C * K; ++n) {
      coupled.tendency.momentum[n] =
          coupled.tendency.momentum[n] + physics.horizontal_momentum_mass_kg_m_s2[n];
      coupled.tendency.potential_temperature_mass[n] +=
          physics.potential_temperature_mass_k_kg_m2_s[n];
    }
    physics_diagnostics = physics.diagnostics;
  } else if (config_.physics.kind == PhysicsKind::kSurfaceEnergyBalance) {
    const auto orbit_state =
        evaluate_orbit(*config_.orbit, config_.planet.rotation_rate_rad_s,
                       s.time_s - config_.run.start_time_s);
    surface_energy_tendency(grid_, *surface_boundary_, s.surface_temperature_k, d,
                            s.surface_pressure_pa, config_.planet, *config_.surface,
                            orbit_state, workspace.surface_physics);
    const auto& physics = workspace.surface_physics;
    if (components != nullptr) {
      components->physics.tendency.momentum = physics.horizontal_momentum_mass_kg_m_s2;
      components->physics.tendency.potential_temperature_mass =
          physics.potential_temperature_mass_k_kg_m2_s;
      components->physics.surface_temperature_k_s = physics.surface_temperature_k_s;
    }
    for (std::size_t n = 0; n < C * K; ++n) {
      coupled.tendency.momentum[n] =
          coupled.tendency.momentum[n] + physics.horizontal_momentum_mass_kg_m_s2[n];
      coupled.tendency.potential_temperature_mass[n] +=
          physics.potential_temperature_mass_k_kg_m2_s[n];
    }
    result.surface_temperature_k_s = physics.surface_temperature_k_s;
    surface_diagnostics = physics.diagnostics;
    surface_dt = physics.stable_time_step_s;
  } else if (config_.physics.kind == PhysicsKind::kPlanetaryNewtonian) {
    std::optional<OrbitState> orbit_state;
    if (config_.physics.geometry == ForcingGeometry::kSubstellar)
      orbit_state = evaluate_orbit(*config_.orbit, config_.planet.rotation_rate_rad_s,
                                   s.time_s - config_.run.start_time_s);
    planetary_newtonian_tendency(
        grid_, coordinate_, d, s.surface_pressure_pa, config_.planet,
        config_.physics.geometry, orbit_state ? &*orbit_state : nullptr,
        workspace.atmospheric_physics, workspace.atmospheric_physics_workspace);
    const auto& physics = workspace.atmospheric_physics;
    if (components != nullptr) {
      components->physics.tendency.momentum = physics.horizontal_momentum_mass_kg_m_s2;
      components->physics.tendency.potential_temperature_mass =
          physics.potential_temperature_mass_k_kg_m2_s;
    }
    for (std::size_t n = 0; n < C * K; ++n) {
      coupled.tendency.momentum[n] =
          coupled.tendency.momentum[n] + physics.horizontal_momentum_mass_kg_m_s2[n];
      coupled.tendency.potential_temperature_mass[n] +=
          physics.potential_temperature_mass_k_kg_m2_s[n];
    }
    physics_diagnostics = physics.diagnostics;
  } else if (config_.physics.kind == PhysicsKind::kGrayRadiation &&
             config_.physics_schedule.kind == PhysicsScheduleKind::kLegacy) {
    const auto orbit_state =
        evaluate_orbit(*config_.orbit, config_.planet.rotation_rate_rad_s,
                       s.time_s - config_.run.start_time_s);
    const auto radiation_wall_start = std::chrono::steady_clock::now();
    gray_radiation_tendency(
        grid_, coordinate_, *surface_boundary_, s.surface_temperature_k, d,
        s.surface_pressure_pa, config_.planet, *config_.surface, *config_.radiation,
        orbit_state, config_.boundary_layer.kind == BoundaryLayerKind::kNone,
        workspace.gray_radiation_physics, workspace.gray_radiation_workspace);
    radiation_wall_seconds =
        std::chrono::duration<Real>(std::chrono::steady_clock::now() -
                                    radiation_wall_start)
            .count();
    radiation_column_call_count = C;
    const auto& physics = workspace.gray_radiation_physics;
    if (components != nullptr) {
      components->physics.tendency.momentum = physics.horizontal_momentum_mass_kg_m_s2;
      components->physics.tendency.potential_temperature_mass =
          physics.potential_temperature_mass_k_kg_m2_s;
      components->physics.surface_temperature_k_s = physics.surface_temperature_k_s;
    }
    for (std::size_t n = 0; n < C * K; ++n) {
      coupled.tendency.momentum[n] =
          coupled.tendency.momentum[n] + physics.horizontal_momentum_mass_kg_m_s2[n];
      coupled.tendency.potential_temperature_mass[n] +=
          physics.potential_temperature_mass_k_kg_m2_s[n];
    }
    result.surface_temperature_k_s = physics.surface_temperature_k_s;
    radiation_diagnostics = physics.diagnostics;
    radiation_dt = physics.stable_time_step_s;
    physics_diagnostics.thermal_energy_rate_w =
        physics.diagnostics.dry_thermal_energy_rate_w;
    physics_diagnostics.rayleigh_drag_work_w = physics.diagnostics.rayleigh_drag_work_w;
  }
  finish_region(6);
  result.surface_pressure_pa_s = coupled.surface_pressure_pa_s;
  result.tendency = coupled.tendency;
  result.horizontal_fast_wave_stable_time_step_s = fast_wave_dt;
  result.horizontal_advective_stable_time_step_s = advective_dt;
  result.diffusion_stable_time_step_s = diffusion_dt;
  result.horizontal_stable_time_step_s = std::min(fast_wave_dt, diffusion_dt);
  result.vertical_stable_time_step_s = vertical_dt;
  result.surface_stable_time_step_s = surface_dt;
  result.radiation_stable_time_step_s = radiation_dt;
  result.maximum_continuity_residual_pa_s = coupled.maximum_continuity_residual_pa_s;
  result.physics_diagnostics = physics_diagnostics;
  result.surface_diagnostics = surface_diagnostics;
  result.radiation_diagnostics = radiation_diagnostics;
  result.radiation_column_call_count = radiation_column_call_count;
  result.radiation_wall_seconds = radiation_wall_seconds;
  result.diffusion_kinetic_energy_rate_w = diffusion_rate;
  finish_region(7);
}
void DryHydrostaticDriver::advance(DryHydrostaticState& s, const Real end,
                                   const DryHydrostaticObserver& obs,
                                   const DryHydrostaticCancel& cancel) const {
  const std::uint64_t initial_observer_step = s.step;
  std::uint64_t last_sampled_step = initial_observer_step;
  const auto observe = [&](const DryHydrostaticState& state,
                           const DryHydrostaticStepDiagnostics& step,
                           const bool force_sample = false) {
    if (!obs) return;
    const bool needs_sample = force_sample || state.step == initial_observer_step ||
                              state.step % config_.diagnostics.interval_steps == 0 ||
                              state.time_s >= end;
    if (!needs_sample) {
      obs(state, nullptr, step);
      return;
    }
    diagnose_dry_hydrostatic_state(state, coordinate_, config_.planet,
                                   orography_.surface_geopotential_m2_s2(),
                                   workspace_.derived, workspace_.diagnosis_workspace);
    const auto& derived = workspace_.derived;
    last_sampled_step = state.step;
    auto sampled = step;
    if (config_.physics.kind == PhysicsKind::kHeldSuarez) {
      held_suarez_tendency(grid_, coordinate_, derived, state.surface_pressure_pa,
                           config_.planet, workspace_.atmospheric_physics,
                           workspace_.atmospheric_physics_workspace);
      sampled.physics_rates = workspace_.atmospheric_physics.diagnostics;
    } else if (config_.physics.kind == PhysicsKind::kPlanetaryNewtonian) {
      std::optional<OrbitState> orbit_state;
      if (config_.physics.geometry == ForcingGeometry::kSubstellar)
        orbit_state = evaluate_orbit(*config_.orbit, config_.planet.rotation_rate_rad_s,
                                     state.time_s - config_.run.start_time_s);
      planetary_newtonian_tendency(
          grid_, coordinate_, derived, state.surface_pressure_pa, config_.planet,
          config_.physics.geometry, orbit_state ? &*orbit_state : nullptr,
          workspace_.atmospheric_physics, workspace_.atmospheric_physics_workspace);
      sampled.physics_rates = workspace_.atmospheric_physics.diagnostics;
    } else if (config_.physics.kind == PhysicsKind::kSurfaceEnergyBalance) {
      const auto orbit_state =
          evaluate_orbit(*config_.orbit, config_.planet.rotation_rate_rad_s,
                         state.time_s - config_.run.start_time_s);
      surface_energy_tendency(grid_, *surface_boundary_, state.surface_temperature_k,
                              derived, state.surface_pressure_pa, config_.planet,
                              *config_.surface, orbit_state,
                              workspace_.surface_physics);
      sampled.surface_rates = workspace_.surface_physics.diagnostics;
    } else if (config_.physics.kind == PhysicsKind::kGrayRadiation) {
      const auto orbit_state =
          evaluate_orbit(*config_.orbit, config_.planet.rotation_rate_rad_s,
                         state.time_s - config_.run.start_time_s);
      gray_radiation_tendency(
          grid_, coordinate_, *surface_boundary_, state.surface_temperature_k, derived,
          state.surface_pressure_pa, config_.planet, *config_.surface,
          *config_.radiation, orbit_state,
          config_.boundary_layer.kind == BoundaryLayerKind::kNone,
          workspace_.gray_radiation_physics, workspace_.gray_radiation_workspace);
      sampled.radiation_rates = workspace_.gray_radiation_physics.diagnostics;
      if (step.accepted_time_step_s == 0.0)
        sampled.radiation_stable_time_step_s =
            workspace_.gray_radiation_physics.stable_time_step_s;
      sampled.physics_rates.thermal_energy_rate_w =
          sampled.radiation_rates.dry_thermal_energy_rate_w;
      sampled.physics_rates.rayleigh_drag_work_w =
          sampled.radiation_rates.rayleigh_drag_work_w;
    }
    obs(state, &derived, sampled);
  };
  observe(s, {});
  std::vector<Vec3> centres;
  centres.reserve(grid_.cell_count());
  for (const auto& cell : grid_.cells()) centres.push_back(cell.center);
  std::optional<std::size_t> water_vapor_tracer;
  if (config_.moisture.kind == MoistureKind::kDiluteWater)
    water_vapor_tracer = TracerRegistry(config_.tracers).water_vapor_index();
  const auto project_momentum = [&](DryHydrostaticState& stage) {
    for (std::size_t n = 0; n < stage.horizontal_momentum_mass_kg_m_s.size(); ++n) {
      const auto cell = n / coordinate_.levels();
      stage.horizontal_momentum_mass_kg_m_s[n] = project_tangent(
          stage.horizontal_momentum_mass_kg_m_s[n], grid_.cells()[cell].center);
    }
  };
  const auto diagnose_and_validate = [&](const DryHydrostaticState& stage) {
    diagnose_dry_hydrostatic_state(stage, coordinate_, config_.planet,
                                   orography_.surface_geopotential_m2_s2(),
                                   workspace_.derived, workspace_.diagnosis_workspace);
    validate_dry_hydrostatic_state(stage, workspace_.derived, centres,
                                   config_.vertical.minimum_surface_pressure_pa,
                                   config_.vertical.maximum_surface_pressure_pa,
                                   config_.vertical.temperature_floor_k);
  };
  const auto apply_local_physics = [&](DryHydrostaticState& stage,
                                       const Real dynamics_time_step_s,
                                       DryMixingStepDiagnostics&
                                           boundary_layer_diagnostics,
                                       DryConvectionDiagnostics& convection_diagnostics,
                                       MoistPhysicsStepDiagnostics& moist_diagnostics,
                                       ScheduledRadiationStep& scheduled_radiation,
                                       std::size_t& boundary_layer_calls,
                                       std::size_t& convection_calls,
                                       std::size_t& physics_substeps,
                                       std::size_t& physics_retries) {
    std::vector<PhysicsEvent> events;
    if (config_.physics_schedule.kind == PhysicsScheduleKind::kLegacy) {
      std::size_t substeps = 1;
      if (config_.moisture.kind == MoistureKind::kDiluteWater) {
        Real maximum_substep = config_.moisture.maximum_physics_substep_s;
        if (config_.convection.kind == ConvectionKind::kSimpleBettsMiller)
          maximum_substep =
              std::min(maximum_substep, 0.25 * config_.convection.relaxation_time_s);
        substeps =
            static_cast<std::size_t>(std::ceil(dynamics_time_step_s / maximum_substep));
      }
      const Real local_time_step = dynamics_time_step_s / static_cast<Real>(substeps);
      events.reserve(substeps);
      for (std::size_t substep = 0; substep < substeps; ++substep)
        events.push_back(
            {.end_time_offset_s = (substep + 1) * local_time_step,
             .interval_s = local_time_step,
             .radiation_interval_s = 0.0,
             .boundary_layer_interval_s = local_time_step,
             .convection_interval_s = local_time_step,
             .radiation = false,
             .boundary_layer = config_.boundary_layer.kind != BoundaryLayerKind::kNone,
             .convection = config_.convection.kind != ConvectionKind::kNone,
             .saturation_adjustment =
                 config_.moisture.kind == MoistureKind::kDiluteWater});
    } else {
      // Radiation still resides in the dynamics RHS in this staged commit. Its
      // deadline joins this event union when the cached-flux update is moved out.
      events =
          make_physics_events(config_.physics_schedule, dynamics_time_step_s,
                              config_.physics.kind == PhysicsKind::kGrayRadiation,
                              config_.boundary_layer.kind != BoundaryLayerKind::kNone,
                              config_.convection.kind != ConvectionKind::kNone,
                              config_.moisture.kind == MoistureKind::kDiluteWater);
    }
    boundary_layer_diagnostics = {};
    convection_diagnostics = {};
    moist_diagnostics = {};
    scheduled_radiation = {};
    if (config_.physics_schedule.convection_update_mode ==
        ConvectionUpdateMode::kCachedRelaxation)
      reset_sbm_reference_cache(workspace_.moist_workspace, grid_.cell_count());
    bool has_convection_diagnostics = false;
    constexpr std::size_t kMaximumPhysicsSubstepRefinements = 10;
    struct PendingPhysicsInterval {
      Real event_interval_s;
      Real end_time_offset_s;
      Real radiation_interval_s;
      Real boundary_layer_interval_s;
      Real convection_interval_s;
      std::size_t refinement;
      bool radiation;
      bool boundary_layer;
      bool convection;
      bool saturation_adjustment;
    };
    DryHydrostaticState saved_stage;
    std::array<PendingPhysicsInterval, kMaximumPhysicsSubstepRefinements + 1>
        pending_intervals{};
    for (const auto& event : events) {
      std::size_t pending_count = 0;
      pending_intervals[pending_count++] = {
          .event_interval_s = event.interval_s,
          .end_time_offset_s = event.end_time_offset_s,
          .radiation_interval_s = event.radiation_interval_s,
          .boundary_layer_interval_s = event.boundary_layer_interval_s,
          .convection_interval_s = event.convection_interval_s,
          .refinement = 0,
          .radiation = event.radiation,
          .boundary_layer = event.boundary_layer,
          .convection = event.convection,
          .saturation_adjustment = event.saturation_adjustment};
      while (pending_count > 0) {
        const auto interval = pending_intervals[--pending_count];
        saved_stage = stage;
        const auto saved_boundary_layer = boundary_layer_diagnostics;
        const auto saved_convection = convection_diagnostics;
        const auto saved_moisture = moist_diagnostics;
        const auto saved_radiation_rates = scheduled_radiation.rates;
        const auto saved_radiation_stable = scheduled_radiation.stable_time_step_s;
        const auto saved_radiation_surface =
            scheduled_radiation.surface_temperature_after_update;
        const auto saved_sbm_cache = workspace_.moist_workspace.convection_cache;
        const bool saved_has_convection = has_convection_diagnostics;
        ++physics_substeps;
        try {
          diagnose_and_validate(stage);
          if (interval.radiation) {
            const Real evaluation_time_s =
                stage.time_s - dynamics_time_step_s + interval.end_time_offset_s;
            const auto orbit_state =
                evaluate_orbit(*config_.orbit, config_.planet.rotation_rate_rad_s,
                               evaluation_time_s - config_.run.start_time_s);
            const auto wall_start = std::chrono::steady_clock::now();
            gray_radiation_tendency(
                grid_, coordinate_, *surface_boundary_, stage.surface_temperature_k,
                workspace_.derived, stage.surface_pressure_pa, config_.planet,
                *config_.surface, *config_.radiation, orbit_state,
                config_.boundary_layer.kind == BoundaryLayerKind::kNone,
                workspace_.gray_radiation_physics, workspace_.gray_radiation_workspace);
            scheduled_radiation.wall_seconds +=
                std::chrono::duration<Real>(std::chrono::steady_clock::now() -
                                            wall_start)
                    .count();
            scheduled_radiation.column_call_count += grid_.cell_count();
            const auto& tendency = workspace_.gray_radiation_physics;
            scheduled_radiation.stable_time_step_s = std::min(
                scheduled_radiation.stable_time_step_s, tendency.stable_time_step_s);
            if (interval.radiation_interval_s > tendency.stable_time_step_s)
              throw std::runtime_error(
                  "scheduled radiation violates its explicit stability limit");
            for (std::size_t n = 0; n < stage.potential_temperature_mass_k_kg_m2.size();
                 ++n) {
              stage.potential_temperature_mass_k_kg_m2[n] +=
                  interval.radiation_interval_s *
                  tendency.potential_temperature_mass_k_kg_m2_s[n];
              stage.horizontal_momentum_mass_kg_m_s[n] =
                  stage.horizontal_momentum_mass_kg_m_s[n] +
                  interval.radiation_interval_s *
                      tendency.horizontal_momentum_mass_kg_m_s2[n];
            }
            for (std::size_t cell = 0; cell < stage.surface_temperature_k.size();
                 ++cell)
              stage.surface_temperature_k[cell] +=
                  interval.radiation_interval_s *
                  tendency.surface_temperature_k_s[cell];
            scheduled_radiation.rates = weighted_radiation_diagnostics(
                scheduled_radiation.rates, 1.0, tendency.diagnostics,
                interval.radiation_interval_s / dynamics_time_step_s);
            scheduled_radiation.surface_temperature_after_update =
                stage.surface_temperature_k;
            project_momentum(stage);
            diagnose_and_validate(stage);
          }
          if (interval.boundary_layer) {
            std::optional<MoistBoundaryLayerCoupling> moist_boundary;
            if (config_.moisture.kind == MoistureKind::kDiluteWater &&
                config_.moisture.surface_exchange == SurfaceMoistureExchange::kBulk) {
              moist_boundary = MoistBoundaryLayerCoupling{
                  .water_vapor_tracer = *water_vapor_tracer,
                  .land_water_kg_m2 = stage.land_water_kg_m2,
                  .bucket_capacity_kg_m2 = config_.surface->hydrology_capacity_kg_m2,
                  .thermodynamics = {
                      .gas_constant_dry_air_j_kg_k = config_.planet.gas_constant_j_kg_k,
                      .heat_capacity_cp_j_kg_k =
                          config_.planet.heat_capacity_cp_j_kg_k}};
            }
            DryMixingStepDiagnostics local_boundary_layer;
            apply_dry_boundary_layer(
                grid_, coordinate_, *surface_boundary_, stage, workspace_.derived,
                config_.planet, *config_.surface, config_.boundary_layer,
                interval.boundary_layer_interval_s, local_boundary_layer,
                workspace_.dry_mixing_workspace,
                moist_boundary ? &*moist_boundary : nullptr);
            boundary_layer_calls += grid_.cell_count();
            accumulate_boundary_layer_diagnostics(
                boundary_layer_diagnostics, local_boundary_layer,
                interval.boundary_layer_interval_s / dynamics_time_step_s);
            if (moist_boundary) {
              stage.cumulative_evaporation_kg += local_boundary_layer.evaporation_kg;
              stage.cumulative_runoff_kg += local_boundary_layer.runoff_kg;
              stage.cumulative_ocean_water_change_kg +=
                  local_boundary_layer.ocean_water_change_kg;
              stage.cumulative_external_outflow_kg +=
                  local_boundary_layer.external_outflow_kg;
              moist_diagnostics.evaporation_kg += local_boundary_layer.evaporation_kg;
              moist_diagnostics.runoff_kg += local_boundary_layer.runoff_kg;
              moist_diagnostics.ocean_water_change_kg +=
                  local_boundary_layer.ocean_water_change_kg;
              moist_diagnostics.external_outflow_kg +=
                  local_boundary_layer.external_outflow_kg;
              moist_diagnostics.water_budget_residual_kg +=
                  local_boundary_layer.water_budget_residual_kg;
              moist_diagnostics.moist_enthalpy_budget_residual_j +=
                  local_boundary_layer.moist_enthalpy_budget_residual_j;
            }
            diagnose_and_validate(stage);
          }
          DryConvectionDiagnostics local_convection;
          if (interval.boundary_layer || interval.convection) {
            local_convection =
                apply_dry_convection(config_, grid_, stage, workspace_.derived,
                                     workspace_.convective_adjustment,
                                     workspace_.convective_adjustment_workspace);
          }
          if (config_.convection.kind != ConvectionKind::kNone &&
              (interval.boundary_layer || interval.convection)) {
            convection_calls += grid_.cell_count();
            accumulate_convection_diagnostics(
                convection_diagnostics, local_convection,
                interval.event_interval_s / dynamics_time_step_s,
                has_convection_diagnostics);
            diagnose_and_validate(stage);
          }
          if (interval.saturation_adjustment) {
            SbmExecution sbm_execution = SbmExecution::kSkip;
            if (interval.convection)
              sbm_execution = config_.physics_schedule.convection_update_mode ==
                                      ConvectionUpdateMode::kCachedRelaxation
                                  ? SbmExecution::kDiagnoseCacheAndApply
                                  : SbmExecution::kDiagnoseAndApply;
            else if (config_.physics_schedule.convection_update_mode ==
                     ConvectionUpdateMode::kCachedRelaxation)
              sbm_execution = SbmExecution::kApplyCached;
            MoistPhysicsStepDiagnostics local;
            apply_moist_column_physics(
                grid_, coordinate_, *surface_boundary_, stage, workspace_.derived,
                *water_vapor_tracer, config_.planet, config_.moisture,
                config_.convection, *config_.surface,
                interval.convection ? interval.convection_interval_s
                                    : interval.event_interval_s,
                local, workspace_.moist_workspace, sbm_execution);
            accumulate_moist_physics_diagnostics(moist_diagnostics, local);
            diagnose_and_validate(stage);
          }
        } catch (const std::runtime_error&) {
          stage = saved_stage;
          boundary_layer_diagnostics = saved_boundary_layer;
          convection_diagnostics = saved_convection;
          moist_diagnostics = saved_moisture;
          scheduled_radiation.rates = saved_radiation_rates;
          scheduled_radiation.stable_time_step_s = saved_radiation_stable;
          scheduled_radiation.surface_temperature_after_update =
              saved_radiation_surface;
          workspace_.moist_workspace.convection_cache = saved_sbm_cache;
          has_convection_diagnostics = saved_has_convection;
          if (interval.refinement == kMaximumPhysicsSubstepRefinements) throw;
          ++physics_retries;
          const auto half = [&](PendingPhysicsInterval part) {
            part.event_interval_s *= 0.5;
            part.radiation_interval_s *= 0.5;
            part.boundary_layer_interval_s *= 0.5;
            part.convection_interval_s *= 0.5;
            ++part.refinement;
            return part;
          };
          auto first = half(interval);
          auto second = first;
          const Real split_time = interval.end_time_offset_s - first.event_interval_s;
          first.end_time_offset_s = split_time;
          second.end_time_offset_s = interval.end_time_offset_s;
          second.event_interval_s = interval.event_interval_s - first.event_interval_s;
          second.radiation_interval_s =
              interval.radiation_interval_s - first.radiation_interval_s;
          second.boundary_layer_interval_s =
              interval.boundary_layer_interval_s - first.boundary_layer_interval_s;
          second.convection_interval_s =
              interval.convection_interval_s - first.convection_interval_s;
          pending_intervals[pending_count++] = second;
          pending_intervals[pending_count++] = first;
        }
      }
    }
  };
  if (config_.dry_hydrostatic.time_integrator !=
      DryHydrostaticTimeIntegrator::kExplicitSspRk3) {
    const auto& parameters = *config_.semi_implicit;
    const bool use_ark2_comparison = config_.dry_hydrostatic.time_integrator ==
                                     DryHydrostaticTimeIntegrator::kArk2ImexComparison;
    const bool reference_is_per_step =
        parameters.reference_update == SemiImplicitReferenceUpdate::kPerStep;
    const GmresOptions linear_options{
        .restart = static_cast<std::size_t>(parameters.gmres_restart),
        .maximum_iterations =
            static_cast<std::size_t>(parameters.linear_maximum_iterations),
        .relative_tolerance = parameters.linear_relative_tolerance,
        .absolute_tolerance = parameters.linear_absolute_tolerance};
    DryHydrostaticRhs initial_rhs;
    DryHydrostaticRhs candidate_rhs;
    std::array<DryHydrostaticRhs, 2> ark2_stage_rhs;
    bool initial_rhs_is_current = false;
    const bool rhs_is_time_independent =
        config_.physics.kind != PhysicsKind::kPlanetaryNewtonian &&
        config_.physics.kind != PhysicsKind::kSurfaceEnergyBalance &&
        (config_.physics.kind != PhysicsKind::kGrayRadiation ||
         config_.physics_schedule.kind == PhysicsScheduleKind::kProcessIntervals);
    while (s.time_s < end) {
      if (cancel && cancel()) {
        // A cancellation can arrive between observations. Flush the pending
        // diagnostic interval without counting an accepted step a second time.
        if (last_sampled_step != s.step) observe(s, {}, true);
        return;
      }
      const auto step_wall_start = std::chrono::steady_clock::now();
      FullRhsEvaluationCounts full_rhs{};
      Real rhs_wall_seconds = 0.0;
      Real linear_solve_wall_seconds = 0.0;
      std::size_t linear_iterations_total = 0;
      std::size_t linear_iterations_maximum = 0;
      Real linear_relative_residual_maximum = 0.0;
      std::size_t accepted_selected_modes = 0;
      std::size_t split_fast_operator_evaluations = 0;
      std::size_t accepted_nonlinear_iterations = 0;
      Real accepted_nonlinear_residual = 0.0;
      std::size_t retry_count = 0;
      std::size_t cfl_retry_count = 0;
      std::size_t invariant_retry_count = 0;
      std::size_t solver_retry_count = 0;
      std::size_t radiation_column_call_count = 0;
      std::size_t boundary_layer_column_call_count = 0;
      std::size_t convection_column_call_count = 0;
      std::size_t physics_substep_count = 0;
      std::size_t physics_retry_count = 0;
      Real radiation_wall_seconds = 0.0;
      DryMixingStepDiagnostics accepted_boundary_layer{};
      DryConvectionDiagnostics accepted_convection{};
      MoistPhysicsStepDiagnostics accepted_moisture{};
      ScheduledRadiationStep scheduled_radiation{};
      std::vector<Real> accepted_radiation_surface_temperature;
      const DryHydrostaticState initial = s;
      if (reference_is_per_step) update_semi_implicit_reference(initial);
      const auto& reference = *semi_implicit_reference_column_;
      const auto& fast_operator = *semi_implicit_fast_operator_;
      const auto& vertical_modes = *semi_implicit_vertical_modes_;
      const auto& external_operator = *semi_implicit_external_operator_;
      const Real requested_dt = std::min(config_.run.time_step_s, end - initial.time_s);
      const auto timed_rhs = [&](const DryHydrostaticState& state,
                                 DryHydrostaticRhs& result) {
        const auto start = std::chrono::steady_clock::now();
        rhs_with_components(state, result, nullptr, false);
        radiation_column_call_count += result.radiation_column_call_count;
        radiation_wall_seconds += result.radiation_wall_seconds;
        rhs_wall_seconds +=
            std::chrono::duration<Real>(std::chrono::steady_clock::now() - start)
                .count();
      };
      if (!initial_rhs_is_current) {
        ++full_rhs.initial;
        timed_rhs(initial, initial_rhs);
      }
      Real dt = std::min({config_.run.time_step_s,
                          semi_implicit_explicit_limit(config_, initial_rhs),
                          end - initial.time_s});
      dt = pressure_limited_time_step(config_, initial, initial_rhs, dt);
      if (!(dt > 0.0) || !std::isfinite(dt))
        throw std::runtime_error("semi-implicit dry time step is invalid");
      if (dt < parameters.minimum_time_step_s && dt < end - initial.time_s)
        throw std::runtime_error(
            "semi-implicit dry time step is below the configured minimum");
      const Real minimum_retry_time_step_s =
          std::min(parameters.minimum_time_step_s, end - initial.time_s);
      const auto attempt_step = [&](const Real attempted_dt) {
        std::vector<std::size_t> selected_modes;
        if (use_ark2_comparison) {
          // The comparator defines L as the complete reference-linear fast operator.
          // Unlike ICI, it has no outer iteration to converge modes omitted by the
          // approximate inverse, so every mode must be solved.
          selected_modes.resize(vertical_modes.mode_count());
          for (std::size_t mode = 0; mode < selected_modes.size(); ++mode)
            selected_modes[mode] = mode;
        } else {
          selected_modes = select_implicit_vertical_modes(
              grid_, vertical_modes, attempted_dt, parameters.wave_cfl_threshold,
              static_cast<std::size_t>(parameters.maximum_implicit_modes));
        }
        accepted_selected_modes = selected_modes.size();
        DryHydrostaticState candidate = initial;
        candidate.time_s = initial.time_s + attempted_dt;
        // ADR 0016 / Benard (2003): an ICI scheme performs a fixed number of
        // quasi-Newton iterations. The residual is a diagnostic, never an acceptance
        // test; only a non-finite or growing residual rejects the attempt.
        Real initial_nonlinear_residual = 0.0;
        const auto record_solve = [&](const DryHydrostaticModalSolveResult& solve) {
          linear_iterations_total += solve.linear_iterations_total;
          linear_iterations_maximum =
              std::max(linear_iterations_maximum, solve.linear_iterations_maximum);
          linear_relative_residual_maximum = std::max(
              linear_relative_residual_maximum, solve.linear_relative_residual_maximum);
          if (!solve.all_converged || solve.equation_residual_norm >
                                          10.0 * parameters.linear_relative_tolerance)
            throw std::runtime_error(
                modal_failure_message(solve.equation_residual_norm));
        };
        if (use_ark2_comparison) {
          accepted_nonlinear_iterations = 0;
          accepted_nonlinear_residual = 0.0;
          const auto table = ark2_imex_table();
          std::array<DryHydrostaticFastTendency, 3> fast_tendencies;
          std::array<const DryHydrostaticRhs*, 3> full_tendencies{&initial_rhs, nullptr,
                                                                  nullptr};
          std::array<const DryHydrostaticFastTendency*, 3> fast_tendency_pointers{
              nullptr, nullptr, nullptr};
          const auto evaluate_fast = [&](const DryHydrostaticState& state,
                                         const std::size_t stage) {
            const auto perturbation =
                make_dry_hydrostatic_fast_perturbation(state, reference);
            apply_dry_hydrostatic_fast_operator(grid_, config_.planet, fast_operator,
                                                perturbation, fast_tendencies[stage],
                                                semi_implicit_workspace_.fast_operator);
            fast_tendency_pointers[stage] = &fast_tendencies[stage];
            ++split_fast_operator_evaluations;
          };
          evaluate_fast(initial, 0);
          for (std::size_t stage = 1; stage < 3; ++stage) {
            candidate = ark2_stage_right_hand_side(
                initial, full_tendencies, fast_tendency_pointers, stage, attempted_dt);
            const auto correction_rhs =
                make_dry_hydrostatic_fast_perturbation(candidate, reference);
            const auto solve_start = std::chrono::steady_clock::now();
            const auto solve = solve_dry_hydrostatic_modal_correction(
                grid_, config_.planet, fast_operator, vertical_modes, selected_modes,
                table.implicit_matrix[stage][stage] * attempted_dt, correction_rhs,
                linear_options, semi_implicit_workspace_.correction,
                semi_implicit_workspace_);
            linear_solve_wall_seconds +=
                std::chrono::duration<Real>(std::chrono::steady_clock::now() -
                                            solve_start)
                    .count();
            record_solve(solve);
            assign_fast_solution(candidate, reference,
                                 semi_implicit_workspace_.correction);
            candidate.time_s = initial.time_s + table.nodes[stage] * attempted_dt;
            project_momentum(candidate);
            if (!pressure_is_in_range(config_, candidate) ||
                !surface_temperature_is_positive(candidate))
              throw std::runtime_error("ARK2 stage violates a prognostic invariant");
            diagnose_and_validate(candidate);
            ++full_rhs.iteration;
            timed_rhs(candidate, ark2_stage_rhs[stage - 1]);
            if (attempted_dt >
                semi_implicit_explicit_limit(config_, ark2_stage_rhs[stage - 1]))
              throw std::runtime_error(
                  "ARK2 stage violates an explicit CFL constraint");
            full_tendencies[stage] = &ark2_stage_rhs[stage - 1];
            evaluate_fast(candidate, stage);
          }
          candidate = ark2_final_update(initial, full_tendencies, attempted_dt);
          candidate.time_s = initial.time_s + attempted_dt;
          project_momentum(candidate);
          if (!pressure_is_in_range(config_, candidate) ||
              !surface_temperature_is_positive(candidate))
            throw std::runtime_error("ARK2 update violates a prognostic invariant");
          diagnose_and_validate(candidate);
        } else {
          accepted_nonlinear_iterations =
              static_cast<std::size_t>(parameters.nonlinear_iterations);
          for (Index iteration = 0; iteration < parameters.nonlinear_iterations;
               ++iteration) {
            const DryHydrostaticRhs* rhs_at_candidate = nullptr;
            if (iteration == 0 && rhs_is_time_independent) {
              rhs_at_candidate = &initial_rhs;
            } else {
              ++full_rhs.iteration;
              timed_rhs(candidate, candidate_rhs);
              rhs_at_candidate = &candidate_rhs;
            }
            const auto residual = crank_nicolson_residual(
                initial, candidate, initial_rhs, *rhs_at_candidate, attempted_dt,
                parameters.implicit_weight, centres, coordinate_.levels());
            const auto scaled_residual = scaled_crank_nicolson_residual(
                residual, initial, reference, external_operator);
            if (!std::isfinite(scaled_residual.norm))
              throw std::runtime_error("Crank-Nicolson residual is not finite");
            if (iteration == 0) initial_nonlinear_residual = scaled_residual.norm;
            accepted_nonlinear_residual = scaled_residual.norm;
            const auto correction_rhs = negative_fast_residual(residual);
            const auto solve_start = std::chrono::steady_clock::now();
            const auto solve = solve_dry_hydrostatic_modal_correction(
                grid_, config_.planet, fast_operator, vertical_modes, selected_modes,
                parameters.implicit_weight * attempted_dt, correction_rhs,
                linear_options, semi_implicit_workspace_.correction,
                semi_implicit_workspace_);
            linear_solve_wall_seconds +=
                std::chrono::duration<Real>(std::chrono::steady_clock::now() -
                                            solve_start)
                    .count();
            record_solve(solve);
            add_semi_implicit_correction(candidate, semi_implicit_workspace_.correction,
                                         residual);
            project_momentum(candidate);
            if (!pressure_is_in_range(config_, candidate) ||
                !surface_temperature_is_positive(candidate))
              throw std::runtime_error("correction violates a prognostic invariant");
            diagnose_and_validate(candidate);
          }
        }
        // Final evaluation at the accepted candidate: it supplies the explicit CFL
        // check, the step energy attribution, the reported residual, and the
        // first-same- as-last tendency reused as the next step's `initial_rhs`.
        ++full_rhs.final;
        timed_rhs(candidate, candidate_rhs);
        if (attempted_dt > semi_implicit_explicit_limit(config_, candidate_rhs))
          throw std::runtime_error("candidate violates an explicit CFL constraint");
        if (!use_ark2_comparison) {
          const auto final_residual = crank_nicolson_residual(
              initial, candidate, initial_rhs, candidate_rhs, attempted_dt,
              parameters.implicit_weight, centres, coordinate_.levels());
          const auto final_scaled = scaled_crank_nicolson_residual(
              final_residual, initial, reference, external_operator);
          accepted_nonlinear_residual = final_scaled.norm;
          if (!std::isfinite(final_scaled.norm) ||
              final_scaled.norm > initial_nonlinear_residual)
            throw std::runtime_error(nonlinear_failure_message(
                initial_nonlinear_residual, final_scaled.norm, final_scaled.rms_norm,
                final_scaled.component, attempted_dt));
        }
        candidate.step = initial.step + 1;
        diagnose_and_validate(candidate);
        if (config_.physics_schedule.kind == PhysicsScheduleKind::kLegacy)
          accepted_radiation_surface_temperature = candidate.surface_temperature_k;
        apply_local_physics(
            candidate, attempted_dt, accepted_boundary_layer, accepted_convection,
            accepted_moisture, scheduled_radiation, boundary_layer_column_call_count,
            convection_column_call_count, physics_substep_count, physics_retry_count);
        if (config_.physics_schedule.kind == PhysicsScheduleKind::kProcessIntervals) {
          radiation_column_call_count += scheduled_radiation.column_call_count;
          radiation_wall_seconds += scheduled_radiation.wall_seconds;
          accepted_radiation_surface_temperature =
              scheduled_radiation.surface_temperature_after_update;
        }
        return candidate;
      };

      DryHydrostaticState candidate;
      while (true) {
        const auto attempt_rhs_begin = full_rhs.total();
        try {
          candidate = attempt_step(dt);
          break;
        } catch (const std::runtime_error& error) {
          full_rhs.discarded += full_rhs.total() - attempt_rhs_begin;
          const std::string_view reason(error.what());
          if (reason.find("CFL") != std::string_view::npos)
            ++cfl_retry_count;
          else if (reason.find("invariant") != std::string_view::npos ||
                   reason.find("pressure") != std::string_view::npos)
            ++invariant_retry_count;
          else
            ++solver_retry_count;
          const Real retry_time_step_s = 0.5 * dt;
          if (!(retry_time_step_s > 0.0) ||
              initial.time_s + retry_time_step_s == initial.time_s ||
              retry_time_step_s < minimum_retry_time_step_s) {
            throw std::runtime_error(
                "semi-implicit dry step failed at the configured minimum: " +
                std::string(error.what()));
          }
          dt = retry_time_step_s;
          ++retry_count;
        }
      }

      s = std::move(candidate);
      const Real explicit_weight = 1.0 - parameters.implicit_weight;
      const auto ark2_weights = ark2_imex_table().explicit_weights;
      const auto ark2_weighted = [&](const auto& value) {
        return ark2_weights[0] * value(initial_rhs) +
               ark2_weights[1] * value(ark2_stage_rhs[0]) +
               ark2_weights[2] * value(ark2_stage_rhs[1]);
      };
      const Real thermal_energy_rate =
          config_.physics_schedule.kind == PhysicsScheduleKind::kProcessIntervals
              ? scheduled_radiation.rates.dry_thermal_energy_rate_w
          : use_ark2_comparison
              ? ark2_weighted([](const DryHydrostaticRhs& rhs) {
                  return rhs.physics_diagnostics.thermal_energy_rate_w;
                })
              : explicit_weight *
                        initial_rhs.physics_diagnostics.thermal_energy_rate_w +
                    parameters.implicit_weight *
                        candidate_rhs.physics_diagnostics.thermal_energy_rate_w;
      const Real drag_energy_rate =
          config_.physics_schedule.kind == PhysicsScheduleKind::kProcessIntervals
              ? scheduled_radiation.rates.rayleigh_drag_work_w
          : use_ark2_comparison
              ? ark2_weighted([](const DryHydrostaticRhs& rhs) {
                  return rhs.physics_diagnostics.rayleigh_drag_work_w;
                })
              : explicit_weight * initial_rhs.physics_diagnostics.rayleigh_drag_work_w +
                    parameters.implicit_weight *
                        candidate_rhs.physics_diagnostics.rayleigh_drag_work_w;
      const Real diffusion_energy_rate =
          use_ark2_comparison
              ? ark2_weighted([](const DryHydrostaticRhs& rhs) {
                  return rhs.diffusion_kinetic_energy_rate_w;
                })
              : explicit_weight * initial_rhs.diffusion_kinetic_energy_rate_w +
                    parameters.implicit_weight *
                        candidate_rhs.diffusion_kinetic_energy_rate_w;
      const auto maximum_courant = [&](const auto stable_time) {
        Real result =
            std::max(courant_from_stable_time_step(dt, stable_time.first,
                                                   stable_time.second(initial_rhs)),
                     courant_from_stable_time_step(dt, stable_time.first,
                                                   stable_time.second(candidate_rhs)));
        if (use_ark2_comparison) {
          result = std::max(
              {result,
               courant_from_stable_time_step(dt, stable_time.first,
                                             stable_time.second(ark2_stage_rhs[0])),
               courant_from_stable_time_step(dt, stable_time.first,
                                             stable_time.second(ark2_stage_rhs[1]))});
        }
        return result;
      };
      Real radiation_stable_time_step_s =
          config_.physics_schedule.kind == PhysicsScheduleKind::kProcessIntervals
              ? scheduled_radiation.stable_time_step_s
              : std::min(initial_rhs.radiation_stable_time_step_s,
                         candidate_rhs.radiation_stable_time_step_s);
      if (use_ark2_comparison)
        radiation_stable_time_step_s =
            std::min({radiation_stable_time_step_s,
                      ark2_stage_rhs[0].radiation_stable_time_step_s,
                      ark2_stage_rhs[1].radiation_stable_time_step_s});
      DryHydrostaticStepDiagnostics step{
          .full_rhs = full_rhs,
          .thermal_energy_contribution_j = dt * thermal_energy_rate,
          .rayleigh_drag_energy_contribution_j = dt * drag_energy_rate,
          .boundary_layer = accepted_boundary_layer,
          .convection = accepted_convection,
          .moisture = accepted_moisture,
          .diffusion_energy_contribution_j = dt * diffusion_energy_rate,
          .requested_time_step_s = requested_dt,
          .accepted_time_step_s = dt,
          .radiation_stable_time_step_s = radiation_stable_time_step_s,
          .advective_cfl = maximum_courant(
              std::pair{config_.dry_hydrostatic.cfl,
                        [](const DryHydrostaticRhs& rhs) {
                          return rhs.horizontal_advective_stable_time_step_s;
                        }}),
          .implicit_wave_courant =
              maximum_vertical_mode_courant(grid_, vertical_modes, dt),
          .vertical_cfl =
              maximum_courant(std::pair{config_.vertical.cfl,
                                        [](const DryHydrostaticRhs& rhs) {
                                          return rhs.vertical_stable_time_step_s;
                                        }}),
          .selected_implicit_modes = accepted_selected_modes,
          .split_fast_operator_evaluations = split_fast_operator_evaluations,
          .linear_iterations_total = linear_iterations_total,
          .linear_iterations_maximum = linear_iterations_maximum,
          .linear_relative_residual_maximum = linear_relative_residual_maximum,
          .nonlinear_iterations = accepted_nonlinear_iterations,
          .nonlinear_relative_residual = accepted_nonlinear_residual,
          .retry_count = retry_count,
          .cfl_retry_count = cfl_retry_count,
          .invariant_retry_count = invariant_retry_count,
          .solver_retry_count = solver_retry_count,
          .radiation_column_call_count = radiation_column_call_count,
          .boundary_layer_column_call_count = boundary_layer_column_call_count,
          .convection_column_call_count = convection_column_call_count,
          .physics_substep_count = physics_substep_count,
          .physics_retry_count = physics_retry_count,
          .radiation_wall_seconds = radiation_wall_seconds,
          .wall_seconds_rhs = rhs_wall_seconds,
          .wall_seconds_linear_solve = linear_solve_wall_seconds};
      if (config_.physics.kind == PhysicsKind::kGrayRadiation) {
        if (config_.physics_schedule.kind == PhysicsScheduleKind::kProcessIntervals)
          step.radiation_rates = scheduled_radiation.rates;
        else if (use_ark2_comparison)
          step.radiation_rates = weighted_radiation_diagnostics(
              initial_rhs.radiation_diagnostics, ark2_weights[0],
              ark2_stage_rhs[0].radiation_diagnostics, ark2_weights[1],
              ark2_stage_rhs[1].radiation_diagnostics, ark2_weights[2]);
        else
          step.radiation_rates = weighted_radiation_diagnostics(
              initial_rhs.radiation_diagnostics, explicit_weight,
              candidate_rhs.radiation_diagnostics, parameters.implicit_weight);
      }
      if (config_.physics.kind == PhysicsKind::kGrayRadiation)
        step.radiation_budget = integrate_gray_radiation_budget(
            grid_, *surface_boundary_, *config_.surface, initial.surface_temperature_k,
            accepted_radiation_surface_temperature, dt, step.radiation_rates);
      if (config_.physics.kind == PhysicsKind::kSurfaceEnergyBalance) {
        const auto weighted = [&](const Real SurfaceEnergyDiagnostics::* member) {
          if (use_ark2_comparison)
            return ark2_weighted([&](const DryHydrostaticRhs& rhs) {
              return rhs.surface_diagnostics.*member;
            });
          return explicit_weight * initial_rhs.surface_diagnostics.*member +
                 parameters.implicit_weight * candidate_rhs.surface_diagnostics.*member;
        };
        const SurfaceEnergyDiagnostics rates{
            .absorbed_stellar_power_w =
                weighted(&SurfaceEnergyDiagnostics::absorbed_stellar_power_w),
            .internal_heat_power_w =
                weighted(&SurfaceEnergyDiagnostics::internal_heat_power_w),
            .outgoing_longwave_power_w =
                weighted(&SurfaceEnergyDiagnostics::outgoing_longwave_power_w),
            .sensible_to_atmosphere_power_w =
                weighted(&SurfaceEnergyDiagnostics::sensible_to_atmosphere_power_w),
            .surface_storage_rate_w =
                weighted(&SurfaceEnergyDiagnostics::surface_storage_rate_w)};
        step.surface_budget = integrate_surface_energy_budget(
            grid_, *surface_boundary_, *config_.surface, initial.surface_temperature_k,
            s.surface_temperature_k, dt, rates, rates, rates);
      }
      if (config_.convection.kind == ConvectionKind::kNone &&
          config_.boundary_layer.kind == BoundaryLayerKind::kNone &&
          config_.moisture.kind == MoistureKind::kNone) {
        std::swap(initial_rhs, candidate_rhs);
        initial_rhs_is_current = true;
      } else {
        initial_rhs_is_current = false;
      }
      step.wall_seconds_total = std::chrono::duration<Real>(
                                    std::chrono::steady_clock::now() - step_wall_start)
                                    .count();
      const bool cancelled = cancel && cancel();
      observe(s, step, cancelled);
      if (cancelled) return;
    }
    return;
  }
  DryHydrostaticRhs rhs1;
  DryHydrostaticRhs rhs2;
  DryHydrostaticRhs rhs3;
  while (s.time_s < end) {
    if (cancel && cancel()) {
      // A cancellation can arrive between observations. Flush the pending
      // diagnostic interval without counting an accepted step a second time.
      if (last_sampled_step != s.step) observe(s, {}, true);
      return;
    }
    const DryHydrostaticState initial = s;
    const Real requested_dt = std::min(config_.run.time_step_s, end - initial.time_s);
    std::size_t cfl_retry_count = 0;
    std::size_t invariant_retry_count = 0;
    std::size_t boundary_layer_column_call_count = 0;
    std::size_t convection_column_call_count = 0;
    std::size_t physics_substep_count = 0;
    std::size_t physics_retry_count = 0;
    DryMixingStepDiagnostics accepted_boundary_layer{};
    DryConvectionDiagnostics accepted_convection{};
    MoistPhysicsStepDiagnostics accepted_moisture{};
    ScheduledRadiationStep scheduled_radiation{};
    std::vector<Real> accepted_radiation_surface_temperature;
    FullRhsEvaluationCounts full_rhs{.initial = 1};
    rhs(initial, rhs1);
    std::size_t radiation_column_call_count = rhs1.radiation_column_call_count;
    Real radiation_wall_seconds = rhs1.radiation_wall_seconds;
    Real dt = std::min(
        {config_.run.time_step_s, stable_time_step(rhs1), end - initial.time_s});
    dt = pressure_limited_time_step(config_, initial, rhs1, dt);
    if (!(dt > 0.0) || !std::isfinite(dt)) {
      throw std::runtime_error("dry driver stable time step is invalid");
    }

    while (true) {
      // Charge the previous failed attempt before beginning the next one.
      full_rhs.discarded = full_rhs.iteration;
      auto stage1 = euler_update(initial, rhs1, dt);
      stage1.time_s = initial.time_s + dt;
      if (!pressure_is_in_range(config_, stage1) ||
          !surface_temperature_is_positive(stage1)) {
        ++invariant_retry_count;
        halve_time_step(initial, dt);
        continue;
      }
      project_momentum(stage1);
      ++full_rhs.iteration;
      rhs(stage1, rhs2);
      radiation_column_call_count += rhs2.radiation_column_call_count;
      radiation_wall_seconds += rhs2.radiation_wall_seconds;
      validate_dry_hydrostatic_state(stage1, workspace_.derived, centres,
                                     config_.vertical.minimum_surface_pressure_pa,
                                     config_.vertical.maximum_surface_pressure_pa,
                                     config_.vertical.temperature_floor_k);
      if (dt > stable_time_step(rhs2)) {
        ++cfl_retry_count;
        halve_time_step(initial, dt);
        continue;
      }

      auto stage2 = convex_update(initial, 0.75, stage1, rhs2, 0.25, dt);
      stage2.time_s = initial.time_s + 0.5 * dt;
      if (!pressure_is_in_range(config_, stage2) ||
          !surface_temperature_is_positive(stage2)) {
        ++invariant_retry_count;
        halve_time_step(initial, dt);
        continue;
      }
      project_momentum(stage2);
      ++full_rhs.iteration;
      rhs(stage2, rhs3);
      radiation_column_call_count += rhs3.radiation_column_call_count;
      radiation_wall_seconds += rhs3.radiation_wall_seconds;
      validate_dry_hydrostatic_state(stage2, workspace_.derived, centres,
                                     config_.vertical.minimum_surface_pressure_pa,
                                     config_.vertical.maximum_surface_pressure_pa,
                                     config_.vertical.temperature_floor_k);
      if (dt > stable_time_step(rhs3)) {
        ++cfl_retry_count;
        halve_time_step(initial, dt);
        continue;
      }

      auto next = convex_update(initial, 1.0 / 3.0, stage2, rhs3, 2.0 / 3.0, dt);
      next.time_s = initial.time_s + dt;
      next.step = initial.step + 1;
      if (!pressure_is_in_range(config_, next) ||
          !surface_temperature_is_positive(next)) {
        ++invariant_retry_count;
        halve_time_step(initial, dt);
        continue;
      }
      project_momentum(next);
      try {
        diagnose_and_validate(next);
        if (config_.physics_schedule.kind == PhysicsScheduleKind::kLegacy)
          accepted_radiation_surface_temperature = next.surface_temperature_k;
        apply_local_physics(
            next, dt, accepted_boundary_layer, accepted_convection, accepted_moisture,
            scheduled_radiation, boundary_layer_column_call_count,
            convection_column_call_count, physics_substep_count, physics_retry_count);
        if (config_.physics_schedule.kind == PhysicsScheduleKind::kProcessIntervals) {
          radiation_column_call_count += scheduled_radiation.column_call_count;
          radiation_wall_seconds += scheduled_radiation.wall_seconds;
          accepted_radiation_surface_temperature =
              scheduled_radiation.surface_temperature_after_update;
        }
      } catch (const std::runtime_error&) {
        ++invariant_retry_count;
        halve_time_step(initial, dt);
        continue;
      }
      s = std::move(next);
      DryHydrostaticStepDiagnostics step{
          .full_rhs = full_rhs,
          .thermal_energy_contribution_j =
              dt *
              (config_.physics_schedule.kind == PhysicsScheduleKind::kProcessIntervals
                   ? scheduled_radiation.rates.dry_thermal_energy_rate_w
                   : rhs1.physics_diagnostics.thermal_energy_rate_w / 6.0 +
                         rhs2.physics_diagnostics.thermal_energy_rate_w / 6.0 +
                         2.0 * rhs3.physics_diagnostics.thermal_energy_rate_w / 3.0),
          .rayleigh_drag_energy_contribution_j =
              dt *
              (config_.physics_schedule.kind == PhysicsScheduleKind::kProcessIntervals
                   ? scheduled_radiation.rates.rayleigh_drag_work_w
                   : rhs1.physics_diagnostics.rayleigh_drag_work_w / 6.0 +
                         rhs2.physics_diagnostics.rayleigh_drag_work_w / 6.0 +
                         2.0 * rhs3.physics_diagnostics.rayleigh_drag_work_w / 3.0),
          .boundary_layer = accepted_boundary_layer,
          .convection = accepted_convection,
          .moisture = accepted_moisture,
          .diffusion_energy_contribution_j =
              dt * (rhs1.diffusion_kinetic_energy_rate_w / 6.0 +
                    rhs2.diffusion_kinetic_energy_rate_w / 6.0 +
                    2.0 * rhs3.diffusion_kinetic_energy_rate_w / 3.0),
          .requested_time_step_s = requested_dt,
          .accepted_time_step_s = dt,
          .radiation_stable_time_step_s = std::min(
              {rhs1.radiation_stable_time_step_s, rhs2.radiation_stable_time_step_s,
               rhs3.radiation_stable_time_step_s,
               scheduled_radiation.stable_time_step_s}),
          .retry_count = cfl_retry_count + invariant_retry_count,
          .cfl_retry_count = cfl_retry_count,
          .invariant_retry_count = invariant_retry_count,
          .radiation_column_call_count = radiation_column_call_count,
          .boundary_layer_column_call_count = boundary_layer_column_call_count,
          .convection_column_call_count = convection_column_call_count,
          .physics_substep_count = physics_substep_count,
          .physics_retry_count = physics_retry_count,
          .radiation_wall_seconds = radiation_wall_seconds};
      if (config_.physics.kind == PhysicsKind::kGrayRadiation)
        step.radiation_rates =
            config_.physics_schedule.kind == PhysicsScheduleKind::kProcessIntervals
                ? scheduled_radiation.rates
                : weighted_radiation_diagnostics(rhs1.radiation_diagnostics, 1.0 / 6.0,
                                                 rhs2.radiation_diagnostics, 1.0 / 6.0,
                                                 rhs3.radiation_diagnostics, 2.0 / 3.0);
      if (config_.physics.kind == PhysicsKind::kGrayRadiation)
        step.radiation_budget = integrate_gray_radiation_budget(
            grid_, *surface_boundary_, *config_.surface, initial.surface_temperature_k,
            accepted_radiation_surface_temperature, dt, step.radiation_rates);
      if (config_.physics.kind == PhysicsKind::kSurfaceEnergyBalance) {
        step.surface_budget = integrate_surface_energy_budget(
            grid_, *surface_boundary_, *config_.surface, initial.surface_temperature_k,
            s.surface_temperature_k, dt, rhs1.surface_diagnostics,
            rhs2.surface_diagnostics, rhs3.surface_diagnostics);
      }
      const bool cancelled = cancel && cancel();
      observe(s, step, cancelled);
      if (cancelled) return;
      break;
    }
  }
}
}  // namespace mps
