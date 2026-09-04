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
  Real maximum_wave_speed_m_s;
};
[[nodiscard]] DryHydrostaticEdgeFlux rusanov_dry_hydrostatic_flux(
    const DryHydrostaticPrimitive& left, const DryHydrostaticPrimitive& right,
    const EdgeTangentBasis& basis, Real gas_constant_j_kg_k,
    Real heat_capacity_cp_j_kg_k);
}  // namespace mps
