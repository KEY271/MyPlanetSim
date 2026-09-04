#pragma once
#include <iosfwd>
#include <span>

#include "myplanetsim/dynamics/dry_hydrostatic_sources.hpp"
#include "myplanetsim/dynamics/dry_hydrostatic_state.hpp"
#include "myplanetsim/grid/cubed_sphere_grid.hpp"
namespace mps {
struct DryHydrostaticDiagnostics {
  Real dry_mass_kg;
  Real potential_temperature_mass_k_kg;
  Real tracer_mass_kg;
  Real total_energy_j;
  Real axial_angular_momentum_kg_m2_s;
  Real minimum_surface_pressure_pa;
  Real maximum_surface_pressure_pa;
  Real minimum_temperature_k;
  Real maximum_temperature_k;
};
struct TerrainDiagnostics {
  Real minimum_surface_height_m;
  Real maximum_surface_height_m;
  Real maximum_surface_slope;
  Real pressure_gradient_l1_m_s2;
  Real pressure_gradient_l2_m_s2;
  Real pressure_gradient_linf_m_s2;
  Real maximum_wind_m_s;
  Real model_pressure_gradient_axial_torque_n_m;
  Real boundary_axial_torque_n_m;
  Real axial_torque_residual_n_m;
  Real terrain_pressure_work_w;
};
[[nodiscard]] DryHydrostaticDiagnostics diagnose_dry_hydrostatic_budgets(
    const CubedSphereGrid&, const DryHydrostaticState&, const DryHydrostaticDerived&,
    const PlanetParameters&);
[[nodiscard]] TerrainDiagnostics diagnose_terrain_budgets(
    const CubedSphereGrid&, const DryHydrostaticState&,
    const DryHydrostaticDerived&,
    const DryHydrostaticSources&, std::span<const Real> surface_geopotential_m2_s2,
    const PlanetParameters&);
[[nodiscard]] Real absolute_pressure_velocity_pa_s(
    Real b_half, Real surface_pressure_tendency_pa_s, Vec3 interface_velocity_m_s,
    Vec3 interface_pressure_gradient_pa_m, Real relative_mass_flux_kg_m2_s,
    Real gravity_m_s2);
void write_terrain_diagnostics(std::ostream&, const TerrainDiagnostics&);
}  // namespace mps
