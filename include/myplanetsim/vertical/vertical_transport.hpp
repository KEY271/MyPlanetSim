#pragma once

#include <span>
#include <vector>

#include "myplanetsim/config/experiment_config.hpp"
#include "myplanetsim/core/types.hpp"

namespace mps {

struct VerticalMassFlux {
  Real surface_pressure_tendency_pa_s = 0.0;
  std::vector<Real> target_air_mass_tendency_kg_m2_s;
  std::vector<Real> interface_flux_kg_m2_s;
  Real continuity_residual_pa_s = 0.0;
};

[[nodiscard]] VerticalMassFlux diagnose_vertical_mass_flux(
    const std::vector<Real>& horizontal_air_mass_tendency_kg_m2_s,
    const std::vector<Real>& b_half, Real gravity_m_s2);
void diagnose_vertical_mass_flux(
    std::span<const Real> horizontal_air_mass_tendency_kg_m2_s,
    std::span<const Real> b_half, Real gravity_m_s2, VerticalMassFlux& result);
[[nodiscard]] Real vertical_stable_time_step(
    std::span<const Real> air_mass_kg_m2, std::span<const Real> interface_flux_kg_m2_s,
    std::span<const Real> horizontal_air_mass_tendency_kg_m2_s, Real cfl,
    Real maximum_time_step_s);
[[nodiscard]] Real vertical_maximum_cfl(
    std::span<const Real> air_mass_kg_m2, std::span<const Real> interface_flux_kg_m2_s,
    std::span<const Real> horizontal_air_mass_tendency_kg_m2_s, Real time_step_s);
[[nodiscard]] std::vector<Real> vertical_scalar_rhs(
    const std::vector<Real>& scalar_mass, const std::vector<Real>& air_mass_kg_m2,
    const std::vector<Real>& horizontal_scalar_tendency,
    const VerticalMassFlux& mass_flux, VerticalTransportScheme scheme,
    VerticalLimiterKind limiter);
void vertical_scalar_rhs(std::span<const Real> scalar_mass,
                         std::span<const Real> air_mass_kg_m2,
                         std::span<const Real> horizontal_scalar_tendency,
                         const VerticalMassFlux& mass_flux,
                         VerticalTransportScheme scheme, VerticalLimiterKind limiter,
                         std::span<const Real> interface_coordinate,
                         std::span<const Real> center_coordinate,
                         std::span<Real> scalar_workspace,
                         std::span<Real> flux_workspace, std::span<Real> result);

}  // namespace mps
