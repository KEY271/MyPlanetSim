#pragma once

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
[[nodiscard]] Real vertical_stable_time_step(
    const std::vector<Real>& air_mass_kg_m2,
    const std::vector<Real>& interface_flux_kg_m2_s,
    const std::vector<Real>& horizontal_air_mass_tendency_kg_m2_s, Real cfl,
    Real maximum_time_step_s);
[[nodiscard]] std::vector<Real> vertical_scalar_rhs(
    const std::vector<Real>& scalar_mass, const std::vector<Real>& air_mass_kg_m2,
    const std::vector<Real>& horizontal_scalar_tendency,
    const VerticalMassFlux& mass_flux, VerticalTransportScheme scheme,
    VerticalLimiterKind limiter);

}  // namespace mps
