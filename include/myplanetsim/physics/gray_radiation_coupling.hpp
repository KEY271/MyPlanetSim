#pragma once

#include <span>
#include <vector>

#include "myplanetsim/core/orbit_parameters.hpp"
#include "myplanetsim/dynamics/dry_hydrostatic_state.hpp"
#include "myplanetsim/dynamics/surface_boundary.hpp"
#include "myplanetsim/grid/cubed_sphere_grid.hpp"
#include "myplanetsim/physics/gray_radiation.hpp"

namespace mps {

struct GrayRadiationDiagnostics {
  Real toa_incoming_shortwave_power_w = 0.0;
  Real toa_reflected_shortwave_power_w = 0.0;
  Real toa_outgoing_longwave_power_w = 0.0;
  Real toa_net_upward_power_w = 0.0;
  Real surface_down_shortwave_power_w = 0.0;
  Real surface_up_shortwave_power_w = 0.0;
  Real surface_down_longwave_power_w = 0.0;
  Real surface_up_longwave_power_w = 0.0;
  Real atmospheric_shortwave_heating_power_w = 0.0;
  Real atmospheric_longwave_heating_power_w = 0.0;
  Real surface_storage_rate_w = 0.0;
  Real sensible_to_atmosphere_power_w = 0.0;
  Real internal_heat_power_w = 0.0;
  Real interface_conservation_residual_w = 0.0;
  Real dry_thermal_energy_rate_w = 0.0;
  Real rayleigh_drag_work_w = 0.0;
};

struct GrayRadiationTendency {
  std::vector<Real> surface_temperature_k_s;
  std::vector<Real> potential_temperature_mass_k_kg_m2_s;
  std::vector<Vec3> horizontal_momentum_mass_kg_m_s2;
  GrayRadiationDiagnostics diagnostics;
  Real stable_time_step_s = 0.0;
};

struct GrayRadiationCouplingWorkspace {
  HybridPressureGeometry vertical_geometry;
  GrayRadiativeColumnTendency column_tendency;
  GrayRadiativeColumnWorkspace column_workspace;
};

void gray_radiation_tendency(
    const CubedSphereGrid& grid, const AtmosphericHybridCoordinate& coordinate,
    const SurfaceBoundary& boundary, std::span<const Real> surface_temperature_k,
    const DryHydrostaticDerived& atmosphere, std::span<const Real> surface_pressure_pa,
    const PlanetParameters& planet, const SurfaceParameters& surface,
    const RadiationParameters& radiation, const OrbitState& orbit_state,
    GrayRadiationTendency& result, GrayRadiationCouplingWorkspace& workspace);

}  // namespace mps
