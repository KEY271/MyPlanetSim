#pragma once
#include "myplanetsim/dynamics/dry_hydrostatic_state.hpp"
#include "myplanetsim/vertical/vertical_transport.hpp"
namespace mps {
struct DryHydrostaticTransportTendency {
  std::vector<Real> air_mass;
  std::vector<Vec3> momentum;
  std::vector<Real> potential_temperature_mass;
  std::vector<Real> tracer_mass;
};
struct DryHydrostaticCoupling {
  std::vector<Real> surface_pressure_pa_s;
  DryHydrostaticTransportTendency tendency;
  std::vector<Real> interface_mass_flux_kg_m2_s;
  Real maximum_continuity_residual_pa_s = 0;
};
struct DryHydrostaticCouplingWorkspace {
  VerticalMassFlux mass_flux;
  std::vector<Real> interface_coordinate;
  std::vector<Real> center_coordinate;
  std::vector<Real> scalar;
  std::vector<Real> flux;
  std::vector<Real> rhs;
  std::vector<Real> component;
  std::vector<Real> horizontal_component;
};
[[nodiscard]] DryHydrostaticCoupling couple_dry_hydrostatic_columns(
    const DryHydrostaticState& state, const DryHydrostaticDerived& derived,
    const DryHydrostaticTransportTendency& horizontal, std::span<const Real> b_half,
    Real gravity_m_s2, VerticalTransportScheme scheme, VerticalLimiterKind limiter);
void couple_dry_hydrostatic_columns(
    const DryHydrostaticState& state, const DryHydrostaticDerived& derived,
    const DryHydrostaticTransportTendency& horizontal, std::span<const Real> b_half,
    Real gravity_m_s2, VerticalTransportScheme scheme, VerticalLimiterKind limiter,
    DryHydrostaticCoupling& result, DryHydrostaticCouplingWorkspace& workspace);
}  // namespace mps
