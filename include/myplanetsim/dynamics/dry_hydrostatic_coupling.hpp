#pragma once
#include <exception>

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
struct DryHydrostaticCouplingColumnWorkspace {
  VerticalMassFlux mass_flux;
  std::vector<Real> interface_coordinate;
  std::vector<Real> center_coordinate;
  std::vector<Real> scalar;
  std::vector<Real> flux;
  std::vector<Real> rhs;
  std::vector<Real> component;
  std::vector<Real> horizontal_component;
};
struct DryHydrostaticCouplingWorkspace {
  std::vector<DryHydrostaticCouplingColumnWorkspace> columns;
  std::vector<Real> continuity_residual_pa_s;
  std::vector<std::exception_ptr> failures;
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
// Dry thermodynamics benefits from a smooth vertical theta profile because nonlinear
// limiter switches degrade its cancellation with hydrostatic pressure work. This
// overload lets the driver keep bounded
// momentum/tracer reconstruction while selecting a separate theta limiter.
void couple_dry_hydrostatic_columns(
    const DryHydrostaticState& state, const DryHydrostaticDerived& derived,
    const DryHydrostaticTransportTendency& horizontal, std::span<const Real> b_half,
    Real gravity_m_s2, VerticalTransportScheme scheme, VerticalLimiterKind limiter,
    VerticalLimiterKind potential_temperature_limiter, DryHydrostaticCoupling& result,
    DryHydrostaticCouplingWorkspace& workspace);
}  // namespace mps
