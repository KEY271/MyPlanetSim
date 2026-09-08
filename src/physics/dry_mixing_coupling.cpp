#include "myplanetsim/physics/dry_mixing_coupling.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace mps {

void apply_dry_boundary_layer(
    const CubedSphereGrid& grid, const AtmosphericHybridCoordinate& coordinate,
    const SurfaceBoundary& boundary, DryHydrostaticState& state,
    const DryHydrostaticDerived& derived, const PlanetParameters& planet,
    const SurfaceParameters& surface, const BoundaryLayerParameters& boundary_layer,
    const Real time_step_s, DryMixingStepDiagnostics& diagnostics,
    DryMixingCouplingWorkspace& workspace, const MoistBoundaryLayerCoupling* moisture) {
  const std::size_t cells = derived.cells;
  const std::size_t levels = derived.levels;
  if (cells != grid.cell_count() || coordinate.levels() != levels || levels == 0 ||
      derived.tracer_count == 0 || state.tracer_count != derived.tracer_count ||
      derived.tracer_mixing_ratio.size() != derived.tracer_count * cells * levels ||
      state.surface_temperature_k.size() != cells ||
      boundary.land_fraction().size() != cells ||
      boundary.surface_geopotential_m2_s2().size() != cells)
    throw std::invalid_argument("dry boundary-layer coupling shape mismatch");
  if (moisture != nullptr && (moisture->water_vapor_tracer >= derived.tracer_count ||
                              moisture->land_water_kg_m2.size() != cells))
    throw std::invalid_argument("moist boundary-layer coupling shape mismatch");
  if (!(time_step_s > 0.0) || !std::isfinite(time_step_s))
    throw std::invalid_argument("dry boundary-layer time step must be positive");

  diagnostics = {};
  Real total_area = 0.0;
  workspace.height_half_m.resize(levels + 1);
  workspace.height_full_m.resize(levels);
  for (std::size_t cell = 0; cell < cells; ++cell) {
    const std::size_t begin = cell * levels;
    const auto theta =
        std::span<const Real>(derived.potential_temperature_k).subspan(begin, levels);
    const auto temperature =
        std::span<const Real>(derived.temperature_k).subspan(begin, levels);
    const auto velocity =
        std::span<const Vec3>(derived.velocity_m_s).subspan(begin, levels);
    coordinate.geometry(state.surface_pressure_pa[cell], planet.gravity_m_s2,
                        planet.gas_constant_j_kg_k, planet.heat_capacity_cp_j_kg_k,
                        planet.reference_pressure_pa, workspace.vertical_geometry);
    const Real surface_geopotential = boundary.surface_geopotential_m2_s2()[cell];
    integrate_hydrostatic_column(workspace.vertical_geometry, theta,
                                 planet.heat_capacity_cp_j_kg_k, planet.gravity_m_s2,
                                 surface_geopotential, workspace.hydrostatic_column);
    const Real surface_height = surface_geopotential / planet.gravity_m_s2;
    for (std::size_t interface = 0; interface <= levels; ++interface)
      workspace.height_half_m[interface] =
          workspace.hydrostatic_column.height_half_m[interface] - surface_height;
    workspace.height_half_m.back() = 0.0;
    for (std::size_t level = 0; level < levels; ++level)
      workspace.height_full_m[level] =
          workspace.hydrostatic_column.height_full_m[level] - surface_height;

    const Real land_fraction = boundary.land_fraction()[cell];
    diagnose_boundary_layer_column(
        {.potential_temperature_k = theta,
         .temperature_k = temperature,
         .velocity_m_s = velocity,
         .pressure_half_pa = workspace.vertical_geometry.pressure_half_pa,
         .exner_half = workspace.vertical_geometry.exner_half,
         .height_half_m = workspace.height_half_m,
         .height_full_m = workspace.height_full_m,
         .surface_temperature_k = state.surface_temperature_k[cell],
         .surface_exner = workspace.vertical_geometry.exner_half.back(),
         .gravity_m_s2 = planet.gravity_m_s2,
         .gas_constant_j_kg_k = planet.gas_constant_j_kg_k,
         .heat_capacity_cp_j_kg_k = planet.heat_capacity_cp_j_kg_k,
         .critical_richardson = boundary_layer.critical_richardson,
         .turbulent_prandtl = boundary_layer.turbulent_prandtl,
         .gustiness_m_s = boundary_layer.gustiness_m_s,
         .land_fraction = land_fraction,
         .land_roughness = {.momentum_m = surface.land_roughness_momentum_m,
                            .heat_m = surface.land_roughness_heat_m},
         .ocean_roughness = {.momentum_m = surface.ocean_roughness_momentum_m,
                             .heat_m = surface.ocean_roughness_heat_m}},
        workspace.bulk);
    const Real surface_capacity =
        mixed_surface_heat_capacity(land_fraction, surface.land_heat_capacity_j_m2_k,
                                    surface.ocean_heat_capacity_j_m2_k);
    Real tracer_mass_change_kg_m2 = 0.0;
    for (std::size_t order = 0; order < derived.tracer_count; ++order) {
      std::size_t tracer_index = order;
      if (moisture != nullptr) {
        if (order == moisture->water_vapor_tracer)
          tracer_index = derived.tracer_count - 1;
        else if (order == derived.tracer_count - 1)
          tracer_index = moisture->water_vapor_tracer;
      }
      const bool is_water =
          moisture != nullptr && tracer_index == moisture->water_vapor_tracer;
      const auto tracer = std::span<const Real>(derived.tracer_mixing_ratio)
                              .subspan(tracer_index * cells * levels + begin, levels);
      implicit_boundary_layer_column(
          {.potential_temperature_k = theta,
           .velocity_m_s = velocity,
           .tracer_mixing_ratio = tracer,
           .air_mass_kg_m2 = workspace.vertical_geometry.air_mass_kg_m2,
           .exner_full = workspace.vertical_geometry.exner_full,
           .exner_half = workspace.vertical_geometry.exner_half,
           .height_full_m = workspace.height_full_m,
           .density_half_kg_m3 = workspace.bulk.density_half_kg_m3,
           .eddy_diffusivity_momentum_m2_s =
               workspace.bulk.eddy_diffusivity_momentum_m2_s,
           .eddy_diffusivity_heat_m2_s = workspace.bulk.eddy_diffusivity_heat_m2_s,
           .eddy_diffusivity_tracer_m2_s = workspace.bulk.eddy_diffusivity_tracer_m2_s,
           .surface_temperature_k = state.surface_temperature_k[cell],
           .surface_exner = workspace.vertical_geometry.exner_half.back(),
           .surface_heat_capacity_j_m2_k = surface_capacity,
           .surface_heat_conductance_w_m2_k =
               workspace.bulk.surface_heat_conductance_w_m2_k,
           .surface_drag_conductance_kg_m2_s =
               workspace.bulk.surface_drag_conductance_kg_m2_s,
           .heat_capacity_cp_j_kg_k = planet.heat_capacity_cp_j_kg_k,
           .time_step_s = time_step_s,
           .enable_surface_water_exchange = is_water,
           .surface_pressure_pa = state.surface_pressure_pa[cell],
           .land_fraction = land_fraction,
           .land_water_kg_m2 = is_water ? moisture->land_water_kg_m2[cell] : 0.0,
           .bucket_capacity_kg_m2 = is_water ? moisture->bucket_capacity_kg_m2 : 0.0,
           .bucket_wet_threshold_fraction = 0.75,
           .surface_water_conductance_land_kg_m2_s =
               is_water ? workspace.bulk.surface_water_conductance_land_kg_m2_s : 0.0,
           .surface_water_conductance_ocean_kg_m2_s =
               is_water ? workspace.bulk.surface_water_conductance_ocean_kg_m2_s : 0.0,
           .moist_thermodynamics =
               is_water ? moisture->thermodynamics : DiluteMoistThermodynamics{}},
          workspace.column, workspace.column_workspace);
      tracer_mass_change_kg_m2 += workspace.column.diagnostics.tracer_mass_change_kg_m2;
      for (std::size_t level = 0; level < levels; ++level) {
        const Real mass = workspace.vertical_geometry.air_mass_kg_m2[level];
        const auto q =
            dry_hydrostatic_tracer_offset(tracer_index, cell, level, cells, levels);
        state.tracer_mass_kg_m2[q] = mass * workspace.column.tracer_mixing_ratio[level];
      }
      if (is_water)
        moisture->land_water_kg_m2[cell] = workspace.column.land_water_kg_m2;
    }

    for (std::size_t level = 0; level < levels; ++level) {
      const Real mass = workspace.vertical_geometry.air_mass_kg_m2[level];
      state.potential_temperature_mass_k_kg_m2[begin + level] =
          mass * workspace.column.potential_temperature_k[level];
      state.horizontal_momentum_mass_kg_m_s[begin + level] =
          mass * workspace.column.velocity_m_s[level];
    }
    state.surface_temperature_k[cell] = workspace.column.surface_temperature_k;

    const Real area = grid.cells()[cell].area_m2;
    total_area += area;
    const auto& bulk = workspace.bulk.diagnostics;
    const auto& column = workspace.column.diagnostics;
    diagnostics.sensible_to_atmosphere_energy_j +=
        area * time_step_s * workspace.column.heat_flux_w_m2.back();
    diagnostics.surface_stress_impulse_magnitude_n_s +=
        area * time_step_s * norm(workspace.column.momentum_flux_kg_m_s2.back());
    diagnostics.mean_boundary_layer_height_m += area * bulk.boundary_layer_height_m;
    diagnostics.maximum_momentum_diffusivity_m2_s =
        std::max(diagnostics.maximum_momentum_diffusivity_m2_s,
                 bulk.maximum_momentum_diffusivity_m2_s);
    diagnostics.maximum_heat_diffusivity_m2_s = std::max(
        diagnostics.maximum_heat_diffusivity_m2_s, bulk.maximum_heat_diffusivity_m2_s);
    if (bulk.shallow_stable_layer_unresolved)
      diagnostics.shallow_unresolved_area_fraction += area;
    if (bulk.reaches_model_top) diagnostics.model_top_area_fraction += area;
    diagnostics.atmospheric_heat_change_j += area * column.atmospheric_heat_change_j_m2;
    diagnostics.surface_heat_change_j += area * column.surface_heat_change_j_m2;
    diagnostics.physical_shear_dissipation_j +=
        area * column.physical_shear_dissipation_j_m2;
    diagnostics.physical_surface_drag_dissipation_j +=
        area * column.physical_surface_drag_dissipation_j_m2;
    diagnostics.backward_euler_dissipation_j +=
        area * column.backward_euler_dissipation_j_m2;
    diagnostics.returned_dissipation_heat_j +=
        area * column.returned_dissipation_heat_j_m2;
    diagnostics.heat_budget_residual_j += area * column.heat_budget_residual_j_m2;
    diagnostics.kinetic_energy_change_j += area * column.kinetic_energy_change_j_m2;
    diagnostics.kinetic_energy_identity_residual_j +=
        area * column.kinetic_energy_identity_residual_j_m2;
    diagnostics.momentum_budget_residual_n_s +=
        area * norm(column.momentum_budget_residual_kg_m_s);
    diagnostics.tracer_mass_change_kg += area * tracer_mass_change_kg_m2;
    diagnostics.evaporation_kg += area * column.evaporation_kg_m2;
    diagnostics.runoff_kg += area * land_fraction * column.runoff_kg_m2_land;
    diagnostics.ocean_water_change_kg += area * column.ocean_water_change_kg_m2;
    diagnostics.external_outflow_kg += area * column.external_outflow_kg_m2;
    diagnostics.moist_enthalpy_budget_residual_j +=
        area * column.moist_enthalpy_budget_residual_j_m2;
    diagnostics.water_budget_residual_kg += area * column.water_budget_residual_kg_m2;
  }
  if (!(total_area > 0.0))
    throw std::runtime_error("dry boundary-layer grid area is invalid");
  diagnostics.mean_boundary_layer_height_m /= total_area;
  diagnostics.shallow_unresolved_area_fraction /= total_area;
  diagnostics.model_top_area_fraction /= total_area;
}

}  // namespace mps
