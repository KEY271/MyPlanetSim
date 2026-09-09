#include "myplanetsim/physics/moist_physics_coupling.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#if defined(MPS_ENABLE_OPENMP)
#include <omp.h>
#endif

#include "myplanetsim/physics/saturation_adjustment.hpp"
#include "myplanetsim/physics/surface_hydrology.hpp"

namespace mps {
namespace {

[[nodiscard]] std::size_t physics_worker_count() {
#if defined(MPS_ENABLE_OPENMP)
  return static_cast<std::size_t>(std::max(1, omp_get_max_threads()));
#else
  return 1;
#endif
}

[[nodiscard]] std::size_t physics_worker_index() {
#if defined(MPS_ENABLE_OPENMP)
  return static_cast<std::size_t>(omp_get_thread_num());
#else
  return 0;
#endif
}

void accumulate_diagnostics(MoistPhysicsStepDiagnostics& total,
                            const MoistPhysicsStepDiagnostics& column) {
  total.evaporation_kg += column.evaporation_kg;
  total.convective_precipitation_kg += column.convective_precipitation_kg;
  total.grid_scale_precipitation_kg += column.grid_scale_precipitation_kg;
  total.runoff_kg += column.runoff_kg;
  total.ocean_water_change_kg += column.ocean_water_change_kg;
  total.external_outflow_kg += column.external_outflow_kg;
  total.water_budget_residual_kg += column.water_budget_residual_kg;
  total.moist_enthalpy_budget_residual_j += column.moist_enthalpy_budget_residual_j;
  total.maximum_relative_humidity =
      std::max(total.maximum_relative_humidity, column.maximum_relative_humidity);
  total.maximum_temperature_increment_k = std::max(
      total.maximum_temperature_increment_k, column.maximum_temperature_increment_k);
  total.maximum_vapor_increment =
      std::max(total.maximum_vapor_increment, column.maximum_vapor_increment);
  total.cape_area_time_integral_j_m2_s_kg += column.cape_area_time_integral_j_m2_s_kg;
  total.cin_area_time_integral_j_m2_s_kg += column.cin_area_time_integral_j_m2_s_kg;
  total.lcl_pressure_area_time_integral_pa_m2_s +=
      column.lcl_pressure_area_time_integral_pa_m2_s;
  total.convection_top_pressure_area_time_integral_pa_m2_s +=
      column.convection_top_pressure_area_time_integral_pa_m2_s;
  total.active_area_time_m2_s += column.active_area_time_m2_s;
  total.deep_area_time_m2_s += column.deep_area_time_m2_s;
  total.shallow_area_time_m2_s += column.shallow_area_time_m2_s;
  total.inactive_area_time_m2_s += column.inactive_area_time_m2_s;
  total.model_top_area_time_m2_s += column.model_top_area_time_m2_s;
  total.deep_column_count += column.deep_column_count;
  total.shallow_column_count += column.shallow_column_count;
  total.inactive_column_count += column.inactive_column_count;
  total.sbm_diagnostic_column_count += column.sbm_diagnostic_column_count;
  total.sbm_relaxation_column_count += column.sbm_relaxation_column_count;
}

}  // namespace

void reset_sbm_reference_cache(MoistPhysicsCouplingWorkspace& workspace,
                               const std::size_t cells) {
  workspace.convection_cache.clear();
  workspace.convection_cache.resize(cells);
}

void apply_moist_column_physics(
    const CubedSphereGrid& grid, const AtmosphericHybridCoordinate& coordinate,
    const SurfaceBoundary& boundary, DryHydrostaticState& state,
    const DryHydrostaticDerived& derived, const std::size_t water_vapor_tracer,
    const PlanetParameters& planet, const MoistureParameters& moisture,
    const ConvectionParameters& convection, const SurfaceParameters& surface,
    const Real time_step_s, MoistPhysicsStepDiagnostics& diagnostics,
    MoistPhysicsCouplingWorkspace& workspace, const SbmExecution sbm_execution) {
  const std::size_t cells = derived.cells;
  const std::size_t levels = derived.levels;
  if (cells != grid.cell_count() || levels == 0 ||
      water_vapor_tracer >= state.tracer_count ||
      state.land_water_kg_m2.size() != cells ||
      boundary.land_fraction().size() != cells || coordinate.levels() != levels ||
      !(time_step_s > 0.0))
    throw std::invalid_argument("moist physics coupling shape mismatch");
  const DiluteMoistThermodynamics thermodynamics{
      .gas_constant_dry_air_j_kg_k = planet.gas_constant_j_kg_k,
      .heat_capacity_cp_j_kg_k = planet.heat_capacity_cp_j_kg_k};
  diagnostics = {};
  const std::size_t workers = physics_worker_count();
  workspace.columns.resize(workers);
  for (auto& column : workspace.columns) column.vapor_mixing_ratio.resize(levels);
  workspace.column_diagnostics.resize(cells);
  workspace.column_failures.resize(cells);
  const auto apply_column = [&](const std::size_t cell,
                                MoistPhysicsStepDiagnostics& diagnostics,
                                MoistPhysicsColumnWorkspace& column_workspace) {
    const std::size_t begin = cell * levels;
    const std::size_t vapor_begin =
        dry_hydrostatic_tracer_offset(water_vapor_tracer, cell, 0, cells, levels);
    coordinate.geometry(state.surface_pressure_pa[cell], planet.gravity_m_s2,
                        planet.gas_constant_j_kg_k, planet.heat_capacity_cp_j_kg_k,
                        planet.reference_pressure_pa,
                        column_workspace.vertical_geometry);
    const auto& mass = column_workspace.vertical_geometry.air_mass_kg_m2;
    const auto& pressure = column_workspace.vertical_geometry.pressure_full_pa;
    const auto& exner = column_workspace.vertical_geometry.exner_full;
    const Real area = grid.cells()[cell].area_m2;
    Real initial_atmospheric_water = 0.0;
    for (std::size_t level = 0; level < levels; ++level)
      initial_atmospheric_water += state.tracer_mass_kg_m2[vapor_begin + level];
    Real convective_rain = 0.0;
    if (convection.kind == ConvectionKind::kSimpleBettsMiller &&
        sbm_execution != SbmExecution::kSkip) {
      for (std::size_t level = 0; level < levels; ++level)
        column_workspace.vapor_mixing_ratio[level] =
            state.tracer_mass_kg_m2[vapor_begin + level] / mass[level];
      const SimpleBettsMillerInput sbm_input{
          .temperature_k =
              std::span<const Real>(derived.temperature_k).subspan(begin, levels),
          .vapor_mixing_ratio = column_workspace.vapor_mixing_ratio,
          .air_mass_kg_m2 = mass,
          .pressure_full_pa = pressure,
          .pressure_half_pa = column_workspace.vertical_geometry.pressure_half_pa,
          .gravity_m_s2 = planet.gravity_m_s2,
          .relative_humidity_reference = convection.reference_relative_humidity,
          .relaxation_time_s = convection.relaxation_time_s,
          .time_step_s = time_step_s,
          .minimum_temperature_k = 150.0,
          .thermodynamics = thermodynamics};
      const SimpleBettsMillerReference* reference = nullptr;
      SbmReferenceCache* cache = nullptr;
      bool reused_reference = false;
      if (sbm_execution == SbmExecution::kDiagnoseAndApply) {
        diagnose_sbm_reference(sbm_input, column_workspace.convection_reference);
        reference = &column_workspace.convection_reference;
        ++diagnostics.sbm_diagnostic_column_count;
      } else {
        if (workspace.convection_cache.size() != cells)
          throw std::invalid_argument("SBM cache is not initialized for this grid");
        cache = &workspace.convection_cache[cell];
        const Real bottom_temperature = sbm_input.temperature_k.back();
        const Real bottom_vapor = sbm_input.vapor_mixing_ratio.back();
        const bool state_invalidates =
            cache->valid &&
            (std::abs(bottom_temperature - cache->bottom_temperature_k) > 0.5 ||
             std::abs(bottom_vapor - cache->bottom_vapor_mixing_ratio) >
                 std::max(0.1 * cache->bottom_vapor_mixing_ratio, 1e-5) ||
             std::abs(state.surface_pressure_pa[cell] - cache->surface_pressure_pa) /
                     cache->surface_pressure_pa >
                 0.01);
        if (sbm_execution == SbmExecution::kDiagnoseCacheAndApply ||
            state_invalidates) {
          diagnose_sbm_reference(sbm_input, cache->reference);
          cache->bottom_temperature_k = bottom_temperature;
          cache->bottom_vapor_mixing_ratio = bottom_vapor;
          cache->surface_pressure_pa = state.surface_pressure_pa[cell];
          cache->valid = true;
          ++diagnostics.sbm_diagnostic_column_count;
        } else {
          reused_reference = cache->valid;
        }
        if (cache->valid) reference = &cache->reference;
      }
      if (reference == nullptr) {
        column_workspace.convection = {};
        column_workspace.convection.temperature_k.assign(
            sbm_input.temperature_k.begin(), sbm_input.temperature_k.end());
        column_workspace.convection.vapor_mixing_ratio.assign(
            sbm_input.vapor_mixing_ratio.begin(), sbm_input.vapor_mixing_ratio.end());
      } else {
        apply_sbm_relaxation(sbm_input, *reference, column_workspace.convection);
        ++diagnostics.sbm_relaxation_column_count;
        // A freshly diagnosed reference balances the column exactly: a deep column
        // dries by the rain it produces and a shallow column conserves its water.
        // A reference kept from an earlier state does not, so the column water
        // residual, not only the bottom-layer change, decides whether the stored
        // reference is still usable.
        if (reused_reference) {
          const auto& applied = column_workspace.convection.diagnostics;
          Real column_water = 0.0;
          for (std::size_t level = 0; level < levels; ++level)
            column_water += state.tracer_mass_kg_m2[vapor_begin + level];
          const Real residual =
              applied.column_water_change_kg_m2 + applied.convective_rain_kg_m2;
          if (std::abs(residual) > 1e-12 * column_water + 1e-12) {
            diagnose_sbm_reference(sbm_input, cache->reference);
            cache->bottom_temperature_k = sbm_input.temperature_k.back();
            cache->bottom_vapor_mixing_ratio = sbm_input.vapor_mixing_ratio.back();
            cache->surface_pressure_pa = state.surface_pressure_pa[cell];
            ++diagnostics.sbm_diagnostic_column_count;
            apply_sbm_relaxation(sbm_input, cache->reference,
                                 column_workspace.convection);
          }
        }
      }
      convective_rain = column_workspace.convection.diagnostics.convective_rain_kg_m2;
      const auto& convective = column_workspace.convection.diagnostics;
      diagnostics.cape_area_time_integral_j_m2_s_kg +=
          area * time_step_s * convective.cape_j_kg;
      diagnostics.cin_area_time_integral_j_m2_s_kg +=
          area * time_step_s * convective.cin_j_kg;
      diagnostics.maximum_temperature_increment_k =
          std::max(diagnostics.maximum_temperature_increment_k,
                   convective.maximum_temperature_increment_k);
      Real maximum_vapor_increment = 0.0;
      for (std::size_t level = 0; level < levels; ++level) {
        const Real initial_vapor =
            state.tracer_mass_kg_m2[vapor_begin + level] / mass[level];
        maximum_vapor_increment =
            std::max(maximum_vapor_increment,
                     std::abs(column_workspace.convection.vapor_mixing_ratio[level] -
                              initial_vapor));
      }
      diagnostics.maximum_vapor_increment =
          std::max(diagnostics.maximum_vapor_increment, maximum_vapor_increment);
      if (convective.branch != MoistConvectionBranch::kNone) {
        diagnostics.active_area_time_m2_s += area * time_step_s;
        diagnostics.lcl_pressure_area_time_integral_pa_m2_s +=
            area * time_step_s * convective.lcl_pressure_pa;
        diagnostics.convection_top_pressure_area_time_integral_pa_m2_s +=
            area * time_step_s * convective.convection_top_pressure_pa;
      }
      if (convective.reaches_model_top)
        diagnostics.model_top_area_time_m2_s += area * time_step_s;
      if (convective.branch == MoistConvectionBranch::kDeep) {
        ++diagnostics.deep_column_count;
        diagnostics.deep_area_time_m2_s += area * time_step_s;
      } else if (convective.branch == MoistConvectionBranch::kShallow) {
        ++diagnostics.shallow_column_count;
        diagnostics.shallow_area_time_m2_s += area * time_step_s;
      } else {
        ++diagnostics.inactive_column_count;
        diagnostics.inactive_area_time_m2_s += area * time_step_s;
      }
      diagnostics.moist_enthalpy_budget_residual_j +=
          area * convective.moist_enthalpy_change_j_m2;
      for (std::size_t level = 0; level < levels; ++level) {
        state.potential_temperature_mass_k_kg_m2[begin + level] =
            mass[level] * column_workspace.convection.temperature_k[level] /
            exner[level];
        state.tracer_mass_kg_m2[vapor_begin + level] =
            mass[level] * column_workspace.convection.vapor_mixing_ratio[level];
      }
    }

    Real grid_scale_rain = 0.0;
    if (moisture.condensation == CondensationKind::kSaturationAdjustment) {
      const auto adjustment = adjust_saturation_column(
          mass, pressure, exner,
          std::span<Real>(state.potential_temperature_mass_k_kg_m2)
              .subspan(begin, levels),
          std::span<Real>(state.tracer_mass_kg_m2).subspan(vapor_begin, levels),
          thermodynamics);
      grid_scale_rain = adjustment.condensed_water_kg_m2;
      diagnostics.moist_enthalpy_budget_residual_j +=
          area * adjustment.enthalpy_change_j_m2;
    }
    const Real precipitation = convective_rain + grid_scale_rain;
    const Real land_fraction = boundary.land_fraction()[cell];
    const auto routed = route_surface_precipitation(land_fraction, precipitation,
                                                    state.land_water_kg_m2[cell],
                                                    surface.hydrology_capacity_kg_m2);
    state.land_water_kg_m2[cell] = routed.final_land_water_kg_m2;
    diagnostics.convective_precipitation_kg += area * convective_rain;
    diagnostics.grid_scale_precipitation_kg += area * grid_scale_rain;
    diagnostics.runoff_kg += area * land_fraction * routed.runoff_kg_m2_land;
    diagnostics.ocean_water_change_kg += area * routed.ocean_gain_kg_m2_cell;
    diagnostics.external_outflow_kg += area * routed.external_outflow_kg_m2_cell;
    Real final_atmospheric_water = 0.0;
    for (std::size_t level = 0; level < levels; ++level) {
      final_atmospheric_water += state.tracer_mass_kg_m2[vapor_begin + level];
      const Real temperature = state.potential_temperature_mass_k_kg_m2[begin + level] /
                               mass[level] * exner[level];
      const Real vapor = state.tracer_mass_kg_m2[vapor_begin + level] / mass[level];
      diagnostics.maximum_relative_humidity = std::max(
          diagnostics.maximum_relative_humidity,
          relative_humidity(vapor, temperature, pressure[level], thermodynamics));
    }
    diagnostics.water_budget_residual_kg +=
        area * ((final_atmospheric_water - initial_atmospheric_water) + precipitation);
  };

#if defined(MPS_ENABLE_OPENMP)
#pragma omp parallel for schedule(static)
#endif
  for (std::ptrdiff_t cell_index = 0; cell_index < static_cast<std::ptrdiff_t>(cells);
       ++cell_index) {
    const auto cell = static_cast<std::size_t>(cell_index);
    workspace.column_diagnostics[cell] = {};
    workspace.column_failures[cell] = nullptr;
    try {
      apply_column(cell, workspace.column_diagnostics[cell],
                   workspace.columns[physics_worker_index()]);
    } catch (...) {
      workspace.column_failures[cell] = std::current_exception();
    }
  }
  for (std::size_t cell = 0; cell < cells; ++cell) {
    if (workspace.column_failures[cell] != nullptr)
      std::rethrow_exception(workspace.column_failures[cell]);
    accumulate_diagnostics(diagnostics, workspace.column_diagnostics[cell]);
  }
  state.cumulative_convective_precipitation_kg +=
      diagnostics.convective_precipitation_kg;
  state.cumulative_grid_scale_precipitation_kg +=
      diagnostics.grid_scale_precipitation_kg;
  state.cumulative_runoff_kg += diagnostics.runoff_kg;
  state.cumulative_ocean_water_change_kg += diagnostics.ocean_water_change_kg;
  state.cumulative_external_outflow_kg += diagnostics.external_outflow_kg;
}

}  // namespace mps
