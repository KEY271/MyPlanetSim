#include "myplanetsim/physics/surface_energy_balance.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "myplanetsim/diagnostics/reductions.hpp"

namespace mps {

SurfaceEnergyTendency surface_energy_tendency(
    const CubedSphereGrid& grid, const SurfaceBoundary& boundary,
    const std::span<const Real> surface_temperature_k,
    const DryHydrostaticDerived& atmosphere,
    const std::span<const Real> surface_pressure_pa, const PlanetParameters& planet,
    const SurfaceParameters& parameters, const OrbitState& orbit_state) {
  SurfaceEnergyTendency result;
  surface_energy_tendency(grid, boundary, surface_temperature_k, atmosphere,
                          surface_pressure_pa, planet, parameters, orbit_state, result);
  return result;
}

void surface_energy_tendency(const CubedSphereGrid& grid,
                             const SurfaceBoundary& boundary,
                             const std::span<const Real> surface_temperature_k,
                             const DryHydrostaticDerived& atmosphere,
                             const std::span<const Real> surface_pressure_pa,
                             const PlanetParameters& planet,
                             const SurfaceParameters& parameters,
                             const OrbitState& orbit_state,
                             SurfaceEnergyTendency& result) {
  const std::size_t volume = atmosphere.cells * atmosphere.levels;
  if (atmosphere.cells != grid.cell_count() || atmosphere.levels == 0 ||
      atmosphere.temperature_k.size() != volume ||
      atmosphere.velocity_m_s.size() != volume ||
      atmosphere.air_mass_kg_m2.size() != volume ||
      surface_pressure_pa.size() != atmosphere.cells ||
      surface_temperature_k.size() != atmosphere.cells ||
      boundary.land_fraction().size() != atmosphere.cells)
    throw std::invalid_argument("surface energy balance shape mismatch");

  result.diagnostics = {};
  result.stable_time_step_s = std::numeric_limits<Real>::infinity();
  result.surface_temperature_k_s.resize(atmosphere.cells);
  result.potential_temperature_mass_k_kg_m2_s.assign(volume, 0.0);
  result.horizontal_momentum_mass_kg_m_s2.assign(volume, {});
  for (std::size_t cell = 0; cell < atmosphere.cells; ++cell) {
    const Real surface_temperature = surface_temperature_k[cell];
    if (!(surface_temperature > 0.0) || !std::isfinite(surface_temperature))
      throw std::invalid_argument("surface temperature must be positive");
    const Real area = grid.cells()[cell].area_m2;
    const Real absorbed = (1.0 - parameters.albedo) * orbit_state.stellar_flux_w_m2 *
                          cosine_solar_zenith(grid.cells()[cell].center, orbit_state);
    const Real surface_temperature_squared = surface_temperature * surface_temperature;
    const Real outgoing = parameters.emissivity * kStefanBoltzmannWm2K4 *
                          surface_temperature_squared * surface_temperature_squared;
    const auto bottom =
        dry_hydrostatic_offset(cell, atmosphere.levels - 1, atmosphere.levels);
    const Real sensible = parameters.air_exchange_coefficient_w_m2_k *
                          (surface_temperature - atmosphere.temperature_k[bottom]);
    const Real capacity = mixed_surface_heat_capacity(
        boundary.land_fraction()[cell], parameters.land_heat_capacity_j_m2_k,
        parameters.ocean_heat_capacity_j_m2_k);
    const Real net =
        absorbed + parameters.internal_heat_flux_w_m2 - outgoing - sensible;
    result.surface_temperature_k_s[cell] = net / capacity;
    const Real exner = std::pow(
        atmosphere.pressure_pa[bottom] / planet.reference_pressure_pa, planet.kappa());
    result.potential_temperature_mass_k_kg_m2_s[bottom] =
        sensible / (planet.heat_capacity_cp_j_kg_k * exner);

    for (std::size_t level = 0; level < atmosphere.levels; ++level) {
      const auto offset = dry_hydrostatic_offset(cell, level, atmosphere.levels);
      const Real sigma = atmosphere.pressure_pa[offset] / surface_pressure_pa[cell];
      const Real drag_rate = std::max(0.0, (sigma - 0.7) / 0.3) / 86400.0;
      result.horizontal_momentum_mass_kg_m_s2[offset] =
          -atmosphere.air_mass_kg_m2[offset] * drag_rate *
          atmosphere.velocity_m_s[offset];
    }

    result.diagnostics.absorbed_stellar_power_w += area * absorbed;
    result.diagnostics.internal_heat_power_w +=
        area * parameters.internal_heat_flux_w_m2;
    result.diagnostics.outgoing_longwave_power_w += area * outgoing;
    result.diagnostics.sensible_to_atmosphere_power_w += area * sensible;
    result.diagnostics.surface_storage_rate_w +=
        area * capacity * result.surface_temperature_k_s[cell];

    const Real longwave_relaxation = 4.0 * parameters.emissivity *
                                     kStefanBoltzmannWm2K4 * surface_temperature *
                                     surface_temperature * surface_temperature;
    if (longwave_relaxation > 0.0)
      result.stable_time_step_s = std::min(
          result.stable_time_step_s, parameters.cfl * capacity / longwave_relaxation);
  }
}

SurfaceEnergyBudget integrate_surface_energy_budget(
    const CubedSphereGrid& grid, const SurfaceBoundary& boundary,
    const SurfaceParameters& parameters,
    const std::span<const Real> initial_surface_temperature_k,
    const std::span<const Real> final_surface_temperature_k, const Real time_step_s,
    const SurfaceEnergyDiagnostics& stage1, const SurfaceEnergyDiagnostics& stage2,
    const SurfaceEnergyDiagnostics& stage3) {
  if (initial_surface_temperature_k.size() != grid.cell_count() ||
      final_surface_temperature_k.size() != grid.cell_count() ||
      boundary.land_fraction().size() != grid.cell_count() || !(time_step_s > 0.0) ||
      !std::isfinite(time_step_s)) {
    throw std::invalid_argument("surface energy budget arguments are invalid");
  }
  const auto integrate = [time_step_s](const Real first, const Real second,
                                       const Real third) {
    return time_step_s * (first / 6.0 + second / 6.0 + 2.0 * third / 3.0);
  };
  SurfaceEnergyBudget budget{
      .absorbed_stellar_energy_j =
          integrate(stage1.absorbed_stellar_power_w, stage2.absorbed_stellar_power_w,
                    stage3.absorbed_stellar_power_w),
      .internal_heat_energy_j =
          integrate(stage1.internal_heat_power_w, stage2.internal_heat_power_w,
                    stage3.internal_heat_power_w),
      .outgoing_longwave_energy_j =
          integrate(stage1.outgoing_longwave_power_w, stage2.outgoing_longwave_power_w,
                    stage3.outgoing_longwave_power_w),
      .sensible_to_atmosphere_energy_j = integrate(
          stage1.sensible_to_atmosphere_power_w, stage2.sensible_to_atmosphere_power_w,
          stage3.sensible_to_atmosphere_power_w)};
  diagnostics::CompensatedAccumulator storage_change;
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
    const Real capacity = mixed_surface_heat_capacity(
        boundary.land_fraction()[cell], parameters.land_heat_capacity_j_m2_k,
        parameters.ocean_heat_capacity_j_m2_k);
    storage_change.add(
        grid.cells()[cell].area_m2 * capacity *
        (final_surface_temperature_k[cell] - initial_surface_temperature_k[cell]));
  }
  budget.surface_storage_change_j = storage_change.value();
  budget.surface_budget_residual_j =
      budget.surface_storage_change_j -
      (budget.absorbed_stellar_energy_j + budget.internal_heat_energy_j -
       budget.outgoing_longwave_energy_j - budget.sensible_to_atmosphere_energy_j);
  return budget;
}

}  // namespace mps
