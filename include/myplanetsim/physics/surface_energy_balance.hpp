#pragma once

#include <span>
#include <vector>

#include "myplanetsim/config/experiment_config.hpp"
#include "myplanetsim/core/orbit_parameters.hpp"
#include "myplanetsim/dynamics/dry_hydrostatic_state.hpp"
#include "myplanetsim/dynamics/surface_boundary.hpp"
#include "myplanetsim/grid/cubed_sphere_grid.hpp"

namespace mps {

inline constexpr Real kStefanBoltzmannWm2K4 = 5.670374419e-8;

struct SurfaceEnergyDiagnostics {
  Real absorbed_stellar_power_w = 0.0;
  Real internal_heat_power_w = 0.0;
  Real outgoing_longwave_power_w = 0.0;
  Real sensible_to_atmosphere_power_w = 0.0;
  Real surface_storage_rate_w = 0.0;
};

struct SurfaceEnergyBudget {
  Real absorbed_stellar_energy_j = 0.0;
  Real internal_heat_energy_j = 0.0;
  Real outgoing_longwave_energy_j = 0.0;
  Real sensible_to_atmosphere_energy_j = 0.0;
  Real surface_storage_change_j = 0.0;
  Real surface_budget_residual_j = 0.0;
};

struct SurfaceEnergyTendency {
  std::vector<Real> surface_temperature_k_s;
  std::vector<Real> potential_temperature_mass_k_kg_m2_s;
  std::vector<Vec3> horizontal_momentum_mass_kg_m_s2;
  SurfaceEnergyDiagnostics diagnostics;
  // cfl * min_c tau_surface(c) with tau = C / (4 eps sigma_SB T_s^3) (ADR 0011). The
  // reservoir is integrated explicitly, so this bounds it before the fact instead of
  // relying on the positivity retry to notice an oscillation.
  Real stable_time_step_s = 0.0;
};

[[nodiscard]] SurfaceEnergyTendency surface_energy_tendency(
    const CubedSphereGrid& grid, const SurfaceBoundary& boundary,
    std::span<const Real> surface_temperature_k,
    const DryHydrostaticDerived& atmosphere, std::span<const Real> surface_pressure_pa,
    const PlanetParameters& planet, const SurfaceParameters& parameters,
    const OrbitState& orbit_state);
void surface_energy_tendency(
    const CubedSphereGrid& grid, const SurfaceBoundary& boundary,
    std::span<const Real> surface_temperature_k,
    const DryHydrostaticDerived& atmosphere, std::span<const Real> surface_pressure_pa,
    const PlanetParameters& planet, const SurfaceParameters& parameters,
    const OrbitState& orbit_state, SurfaceEnergyTendency& result);
[[nodiscard]] SurfaceEnergyBudget integrate_surface_energy_budget(
    const CubedSphereGrid& grid, const SurfaceBoundary& boundary,
    const SurfaceParameters& parameters,
    std::span<const Real> initial_surface_temperature_k,
    std::span<const Real> final_surface_temperature_k, Real time_step_s,
    const SurfaceEnergyDiagnostics& stage1, const SurfaceEnergyDiagnostics& stage2,
    const SurfaceEnergyDiagnostics& stage3);

}  // namespace mps
