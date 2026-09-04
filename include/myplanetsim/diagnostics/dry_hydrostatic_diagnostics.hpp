#pragma once
#include "myplanetsim/dynamics/dry_hydrostatic_state.hpp"
#include "myplanetsim/grid/cubed_sphere_grid.hpp"
namespace mps {
struct DryHydrostaticDiagnostics { Real dry_mass_kg; Real potential_temperature_mass_k_kg; Real tracer_mass_kg; Real total_energy_j; Real axial_angular_momentum_kg_m2_s; Real minimum_surface_pressure_pa; Real maximum_surface_pressure_pa; Real minimum_temperature_k; Real maximum_temperature_k; };
[[nodiscard]] DryHydrostaticDiagnostics diagnose_dry_hydrostatic_budgets(const CubedSphereGrid&,const DryHydrostaticState&,const DryHydrostaticDerived&,const PlanetParameters&);
}
