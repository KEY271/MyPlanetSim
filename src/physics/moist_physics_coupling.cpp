#include "myplanetsim/physics/moist_physics_coupling.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "myplanetsim/physics/saturation_adjustment.hpp"
#include "myplanetsim/physics/surface_hydrology.hpp"

namespace mps {

void apply_moist_column_physics(
    const CubedSphereGrid& grid, const AtmosphericHybridCoordinate& coordinate,
    const SurfaceBoundary& boundary, DryHydrostaticState& state,
    const DryHydrostaticDerived& derived, const std::size_t water_vapor_tracer,
    const PlanetParameters& planet, const MoistureParameters& moisture,
    const ConvectionParameters& convection, const SurfaceParameters& surface,
    const Real time_step_s, MoistPhysicsStepDiagnostics& diagnostics,
    MoistPhysicsCouplingWorkspace& workspace) {
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
  for (std::size_t cell = 0; cell < cells; ++cell) {
    const std::size_t begin = cell * levels;
    const std::size_t vapor_begin =
        dry_hydrostatic_tracer_offset(water_vapor_tracer, cell, 0, cells, levels);
    coordinate.geometry(state.surface_pressure_pa[cell], planet.gravity_m_s2,
                        planet.gas_constant_j_kg_k, planet.heat_capacity_cp_j_kg_k,
                        planet.reference_pressure_pa, workspace.vertical_geometry);
    const auto& mass = workspace.vertical_geometry.air_mass_kg_m2;
    const auto& pressure = workspace.vertical_geometry.pressure_full_pa;
    const auto& exner = workspace.vertical_geometry.exner_full;
    const Real area = grid.cells()[cell].area_m2;
    Real initial_atmospheric_water = 0.0;
    for (std::size_t level = 0; level < levels; ++level)
      initial_atmospheric_water += state.tracer_mass_kg_m2[vapor_begin + level];
    Real convective_rain = 0.0;
    if (convection.kind == ConvectionKind::kSimpleBettsMiller) {
      std::vector<Real> vapor(levels);
      for (std::size_t level = 0; level < levels; ++level)
        vapor[level] = state.tracer_mass_kg_m2[vapor_begin + level] / mass[level];
      const SimpleBettsMillerInput sbm_input{
          .temperature_k =
              std::span<const Real>(derived.temperature_k).subspan(begin, levels),
          .vapor_mixing_ratio = vapor,
          .air_mass_kg_m2 = mass,
          .pressure_full_pa = pressure,
          .pressure_half_pa = workspace.vertical_geometry.pressure_half_pa,
          .gravity_m_s2 = planet.gravity_m_s2,
          .relative_humidity_reference = convection.reference_relative_humidity,
          .relaxation_time_s = convection.relaxation_time_s,
          .time_step_s = time_step_s,
          .minimum_temperature_k = 150.0,
          .thermodynamics = thermodynamics};
      diagnose_sbm_reference(sbm_input, workspace.convection_reference);
      apply_sbm_relaxation(sbm_input, workspace.convection_reference,
                           workspace.convection);
      convective_rain = workspace.convection.diagnostics.convective_rain_kg_m2;
      const auto& convective = workspace.convection.diagnostics;
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
        maximum_vapor_increment = std::max(
            maximum_vapor_increment,
            std::abs(workspace.convection.vapor_mixing_ratio[level] - initial_vapor));
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
            mass[level] * workspace.convection.temperature_k[level] / exner[level];
        state.tracer_mass_kg_m2[vapor_begin + level] =
            mass[level] * workspace.convection.vapor_mixing_ratio[level];
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
