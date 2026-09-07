#include "myplanetsim/physics/gray_radiation_coupling.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "myplanetsim/diagnostics/reductions.hpp"

namespace mps {

void gray_radiation_tendency(
    const CubedSphereGrid& grid, const AtmosphericHybridCoordinate& coordinate,
    const SurfaceBoundary& boundary, const std::span<const Real> surface_temperature_k,
    const DryHydrostaticDerived& atmosphere,
    const std::span<const Real> surface_pressure_pa, const PlanetParameters& planet,
    const SurfaceParameters& surface, const RadiationParameters& radiation,
    const OrbitState& orbit_state, const bool include_legacy_surface_exchange,
    GrayRadiationTendency& result, GrayRadiationCouplingWorkspace& workspace) {
  const std::size_t cells = atmosphere.cells;
  const std::size_t levels = atmosphere.levels;
  const std::size_t volume = cells * levels;
  if (cells != grid.cell_count() || levels == 0 || coordinate.levels() != levels ||
      atmosphere.temperature_k.size() != volume ||
      atmosphere.air_mass_kg_m2.size() != volume ||
      atmosphere.velocity_m_s.size() != volume ||
      atmosphere.exner_full.size() != volume ||
      atmosphere.exner_half.size() != cells * (levels + 1) ||
      surface_pressure_pa.size() != cells || surface_temperature_k.size() != cells ||
      boundary.land_fraction().size() != cells)
    throw std::invalid_argument("gray radiation coupling shape mismatch");

  result.surface_temperature_k_s.resize(cells);
  result.potential_temperature_mass_k_kg_m2_s.assign(volume, 0.0);
  result.horizontal_momentum_mass_kg_m_s2.assign(volume, {});
  result.diagnostics = {};
  result.stable_time_step_s = std::numeric_limits<Real>::infinity();
  diagnostics::CompensatedAccumulator toa_incoming;
  diagnostics::CompensatedAccumulator toa_reflected;
  diagnostics::CompensatedAccumulator toa_outgoing;
  diagnostics::CompensatedAccumulator toa_net;
  diagnostics::CompensatedAccumulator surface_down_sw;
  diagnostics::CompensatedAccumulator surface_up_sw;
  diagnostics::CompensatedAccumulator surface_down_lw;
  diagnostics::CompensatedAccumulator surface_up_lw;
  diagnostics::CompensatedAccumulator atmosphere_sw;
  diagnostics::CompensatedAccumulator atmosphere_lw;
  diagnostics::CompensatedAccumulator surface_storage;
  diagnostics::CompensatedAccumulator sensible_power;
  diagnostics::CompensatedAccumulator internal_power;
  diagnostics::CompensatedAccumulator conservation_residual;
  diagnostics::CompensatedAccumulator thermal_energy_rate;
  diagnostics::CompensatedAccumulator drag_work;

  for (std::size_t cell = 0; cell < cells; ++cell) {
    coordinate.geometry(surface_pressure_pa[cell], planet.gravity_m_s2,
                        planet.gas_constant_j_kg_k, planet.heat_capacity_cp_j_kg_k,
                        planet.reference_pressure_pa, workspace.vertical_geometry);
    const std::size_t begin = cell * levels;
    const std::span<const Real> temperature(atmosphere.temperature_k.data() + begin,
                                            levels);
    const std::span<const Real> exner(atmosphere.exner_full.data() + begin, levels);
    const Real capacity = mixed_surface_heat_capacity(
        boundary.land_fraction()[cell], surface.land_heat_capacity_j_m2_k,
        surface.ocean_heat_capacity_j_m2_k);
    const GrayRadiativeColumnSourceInput input{
        .radiation =
            GrayRadiationColumnInput{
                .pressure_half_pa = workspace.vertical_geometry.pressure_half_pa,
                .temperature_k = temperature,
                .surface_temperature_k = surface_temperature_k[cell],
                .gravity_m_s2 = planet.gravity_m_s2,
                .stellar_flux_w_m2 = orbit_state.stellar_flux_w_m2,
                .cosine_solar_zenith =
                    cosine_solar_zenith(grid.cells()[cell].center, orbit_state),
                .surface_albedo = surface.albedo,
                .surface_emissivity = surface.emissivity,
                .parameters = radiation,
            },
        .exner_full = exner,
        .heat_capacity_cp_j_kg_k = planet.heat_capacity_cp_j_kg_k,
        .surface_heat_capacity_j_m2_k = capacity,
        .air_exchange_coefficient_w_m2_k = include_legacy_surface_exchange
                                               ? surface.air_exchange_coefficient_w_m2_k
                                               : 0.0,
        .internal_heat_flux_w_m2 = surface.internal_heat_flux_w_m2,
    };
    gray_radiative_column_tendency(input, workspace.column_tendency,
                                   workspace.column_workspace);
    const auto& column = workspace.column_tendency;
    result.surface_temperature_k_s[cell] = column.surface_temperature_k_s;
    result.stable_time_step_s =
        std::min(result.stable_time_step_s, column.stable_time_step_s);
    for (std::size_t level = 0; level < levels; ++level) {
      const auto offset = dry_hydrostatic_offset(cell, level, levels);
      result.potential_temperature_mass_k_kg_m2_s[offset] =
          column.potential_temperature_mass_k_kg_m2_s[level];
      if (include_legacy_surface_exchange) {
        const Real sigma = atmosphere.pressure_pa[offset] / surface_pressure_pa[cell];
        const Real drag_rate = std::max(0.0, (sigma - 0.7) / 0.3) / 86400.0;
        result.horizontal_momentum_mass_kg_m_s2[offset] =
            -atmosphere.air_mass_kg_m2[offset] * drag_rate *
            atmosphere.velocity_m_s[offset];
      }
    }

    const Real area = grid.cells()[cell].area_m2;
    const auto& flux = column.fluxes;
    toa_incoming.add(area * flux.shortwave_down_w_m2.front());
    toa_reflected.add(area * flux.shortwave_up_w_m2.front());
    toa_outgoing.add(area * flux.longwave_up_w_m2.front());
    toa_net.add(area * flux.net_flux_w_m2.front());
    surface_down_sw.add(area * flux.shortwave_down_w_m2.back());
    surface_up_sw.add(area * flux.shortwave_up_w_m2.back());
    surface_down_lw.add(area * flux.longwave_down_w_m2.back());
    surface_up_lw.add(area * flux.longwave_up_w_m2.back());
    Real column_sw = 0.0;
    Real column_lw = 0.0;
    Real column_total = 0.0;
    for (std::size_t level = 0; level < levels; ++level) {
      column_sw += flux.shortwave_convergence_w_m2[level];
      column_lw += flux.longwave_convergence_w_m2[level];
      column_total += flux.radiative_convergence_w_m2[level];
    }
    atmosphere_sw.add(area * column_sw);
    atmosphere_lw.add(area * column_lw);
    surface_storage.add(area * column.surface_storage_rate_w_m2);
    sensible_power.add(area * column.sensible_to_atmosphere_w_m2);
    internal_power.add(area * surface.internal_heat_flux_w_m2);
    conservation_residual.add(
        area * (column_total - flux.net_flux_w_m2.back() + flux.net_flux_w_m2.front()));

    Real geopotential_half_rate = 0.0;
    const auto exner_half_begin = cell * (levels + 1);
    for (std::size_t reverse = levels; reverse > 0; --reverse) {
      const std::size_t level = reverse - 1;
      const auto offset = dry_hydrostatic_offset(cell, level, levels);
      const Real theta_rate = column.potential_temperature_mass_k_kg_m2_s[level] /
                              atmosphere.air_mass_kg_m2[offset];
      const Real geopotential_full_rate =
          geopotential_half_rate +
          planet.heat_capacity_cp_j_kg_k * theta_rate *
              (atmosphere.exner_half[exner_half_begin + level + 1] -
               atmosphere.exner_full[offset]);
      thermal_energy_rate.add(
          area * (planet.heat_capacity_cv_j_kg_k() * atmosphere.exner_full[offset] *
                      column.potential_temperature_mass_k_kg_m2_s[level] +
                  atmosphere.air_mass_kg_m2[offset] * geopotential_full_rate));
      geopotential_half_rate += planet.heat_capacity_cp_j_kg_k * theta_rate *
                                (atmosphere.exner_half[exner_half_begin + level + 1] -
                                 atmosphere.exner_half[exner_half_begin + level]);
      drag_work.add(area * dot(atmosphere.velocity_m_s[offset],
                               result.horizontal_momentum_mass_kg_m_s2[offset]));
    }
  }

  result.diagnostics = {
      .toa_incoming_shortwave_power_w = toa_incoming.value(),
      .toa_reflected_shortwave_power_w = toa_reflected.value(),
      .toa_outgoing_longwave_power_w = toa_outgoing.value(),
      .toa_net_upward_power_w = toa_net.value(),
      .surface_down_shortwave_power_w = surface_down_sw.value(),
      .surface_up_shortwave_power_w = surface_up_sw.value(),
      .surface_down_longwave_power_w = surface_down_lw.value(),
      .surface_up_longwave_power_w = surface_up_lw.value(),
      .atmospheric_shortwave_heating_power_w = atmosphere_sw.value(),
      .atmospheric_longwave_heating_power_w = atmosphere_lw.value(),
      .surface_storage_rate_w = surface_storage.value(),
      .sensible_to_atmosphere_power_w = sensible_power.value(),
      .internal_heat_power_w = internal_power.value(),
      .interface_conservation_residual_w = conservation_residual.value(),
      .dry_thermal_energy_rate_w = thermal_energy_rate.value(),
      .rayleigh_drag_work_w = drag_work.value(),
  };
}

GrayRadiationBudget integrate_gray_radiation_budget(
    const CubedSphereGrid& grid, const SurfaceBoundary& boundary,
    const SurfaceParameters& surface,
    const std::span<const Real> initial_surface_temperature_k,
    const std::span<const Real> final_surface_temperature_k, const Real time_step_s,
    const GrayRadiationDiagnostics& rates) {
  if (initial_surface_temperature_k.size() != grid.cell_count() ||
      final_surface_temperature_k.size() != grid.cell_count() ||
      boundary.land_fraction().size() != grid.cell_count() || !(time_step_s > 0.0) ||
      !std::isfinite(time_step_s))
    throw std::invalid_argument("gray radiation budget arguments are invalid");
  diagnostics::CompensatedAccumulator storage_change;
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
    const Real capacity = mixed_surface_heat_capacity(
        boundary.land_fraction()[cell], surface.land_heat_capacity_j_m2_k,
        surface.ocean_heat_capacity_j_m2_k);
    storage_change.add(
        grid.cells()[cell].area_m2 * capacity *
        (final_surface_temperature_k[cell] - initial_surface_temperature_k[cell]));
  }
  GrayRadiationBudget result{
      .toa_incoming_shortwave_energy_j =
          time_step_s * rates.toa_incoming_shortwave_power_w,
      .toa_reflected_shortwave_energy_j =
          time_step_s * rates.toa_reflected_shortwave_power_w,
      .toa_outgoing_longwave_energy_j =
          time_step_s * rates.toa_outgoing_longwave_power_w,
      .toa_net_upward_energy_j = time_step_s * rates.toa_net_upward_power_w,
      .atmospheric_shortwave_heating_energy_j =
          time_step_s * rates.atmospheric_shortwave_heating_power_w,
      .atmospheric_longwave_heating_energy_j =
          time_step_s * rates.atmospheric_longwave_heating_power_w,
      .sensible_to_atmosphere_energy_j =
          time_step_s * rates.sensible_to_atmosphere_power_w,
      .internal_heat_energy_j = time_step_s * rates.internal_heat_power_w,
      .surface_storage_change_j = storage_change.value(),
      .interface_conservation_residual_j =
          time_step_s * rates.interface_conservation_residual_w,
      .dry_thermal_energy_j = time_step_s * rates.dry_thermal_energy_rate_w,
      .rayleigh_drag_energy_j = time_step_s * rates.rayleigh_drag_work_w,
  };
  result.surface_time_integration_residual_j =
      result.surface_storage_change_j - time_step_s * rates.surface_storage_rate_w;
  return result;
}

}  // namespace mps
