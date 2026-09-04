#pragma once
#include "myplanetsim/dynamics/dry_hydrostatic_state.hpp"
#include "myplanetsim/vertical/vertical_transport.hpp"
namespace mps {
struct DryHydrostaticTransportTendency { std::vector<Real> air_mass; std::vector<Vec3> momentum; std::vector<Real> potential_temperature_mass; std::vector<Real> tracer_mass; };
struct DryHydrostaticCoupling { std::vector<Real> surface_pressure_pa_s; DryHydrostaticTransportTendency tendency; std::vector<Real> interface_mass_flux_kg_m2_s; Real maximum_continuity_residual_pa_s=0; };
[[nodiscard]] DryHydrostaticCoupling couple_dry_hydrostatic_columns(const DryHydrostaticState& state,const DryHydrostaticDerived& derived,const DryHydrostaticTransportTendency& horizontal,std::span<const Real> b_half,Real gravity_m_s2,VerticalTransportScheme scheme,VerticalLimiterKind limiter);
}
