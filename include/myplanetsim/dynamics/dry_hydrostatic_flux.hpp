#pragma once
#include "myplanetsim/numerics/spherical_operators.hpp"
namespace mps {
struct DryHydrostaticPrimitive {
  Real air_mass_kg_m2;
  Vec3 velocity_m_s;
  Real potential_temperature_k;
  Real tracer_mixing_ratio;
  Real temperature_k;
};
struct DryHydrostaticEdgeFlux {
  Real air_mass_kg_m_s;
  Vec3 momentum_kg_s2;
  Real potential_temperature_mass_k_kg_m_s;
  Real tracer_mass_kg_m_s;
  // Advective characteristic used only in the jump dissipation. The hydrostatic
  // external/Lamb speed remains in maximum_wave_speed_m_s for the explicit CFL.
  Real maximum_dissipation_speed_m_s;
  Real maximum_wave_speed_m_s;
};
[[nodiscard]] DryHydrostaticEdgeFlux rusanov_dry_hydrostatic_flux(
    const DryHydrostaticPrimitive& left, const DryHydrostaticPrimitive& right,
    const EdgeTangentBasis& basis, Real gas_constant_j_kg_k,
    Real heat_capacity_cp_j_kg_k, bool compute_maximum_wave_speed = true);
}  // namespace mps
