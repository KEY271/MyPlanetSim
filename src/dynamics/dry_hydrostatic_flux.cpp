#include "myplanetsim/dynamics/dry_hydrostatic_flux.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
namespace mps {
DryHydrostaticEdgeFlux rusanov_dry_hydrostatic_flux(const DryHydrostaticPrimitive& l,
                                                    const DryHydrostaticPrimitive& r,
                                                    const EdgeTangentBasis& b,
                                                    const Real rd, const Real cp,
                                                    const bool compute_wave_speed) {
  if (!(l.air_mass_kg_m2 > 0 && r.air_mass_kg_m2 > 0 && l.temperature_k > 0 &&
        r.temperature_k > 0 && cp > rd && rd > 0) ||
      !is_finite(l.velocity_m_s) || !is_finite(r.velocity_m_s))
    throw std::invalid_argument("invalid dry hydrostatic edge state");
  const auto ul = dot(l.velocity_m_s, b.normal) * b.normal +
                  dot(l.velocity_m_s, b.tangent) * b.tangent;
  const auto ur = dot(r.velocity_m_s, b.normal) * b.normal +
                  dot(r.velocity_m_s, b.tangent) * b.tangent;
  const auto unl = dot(ul, b.normal), unr = dot(ur, b.normal);
  // These physical fluxes are the advective block of the split hydrostatic system;
  // pressure is advanced separately by dry_hydrostatic_sources. Upwind their jumps
  // with the advective characteristic, while retaining the fast horizontal Lamb mode
  // below for the SSPRK3 stability limit.
  const auto dissipation_speed = std::max(std::abs(unl), std::abs(unr));
  const auto maximum_wave_speed =
      compute_wave_speed
          ? std::max(std::abs(unl) + std::sqrt(cp / (cp - rd) * rd * l.temperature_k),
                     std::abs(unr) + std::sqrt(cp / (cp - rd) * rd * r.temperature_k))
          : 0.0;
  const auto fm = 0.5 * (l.air_mass_kg_m2 * unl + r.air_mass_kg_m2 * unr) -
                  0.5 * dissipation_speed * (r.air_mass_kg_m2 - l.air_mass_kg_m2);
  const auto ml = l.air_mass_kg_m2 * ul, mr = r.air_mass_kg_m2 * ur;
  return {.air_mass_kg_m_s = fm,
          .momentum_kg_s2 =
              0.5 * (unl * ml + unr * mr) - 0.5 * dissipation_speed * (mr - ml),
          .potential_temperature_mass_k_kg_m_s =
              0.5 * (l.air_mass_kg_m2 * unl * l.potential_temperature_k +
                     r.air_mass_kg_m2 * unr * r.potential_temperature_k) -
              0.5 * dissipation_speed *
                  (r.air_mass_kg_m2 * r.potential_temperature_k -
                   l.air_mass_kg_m2 * l.potential_temperature_k),
          .tracer_mass_kg_m_s = 0.5 * (l.air_mass_kg_m2 * unl * l.tracer_mixing_ratio +
                                       r.air_mass_kg_m2 * unr * r.tracer_mixing_ratio) -
                                0.5 * dissipation_speed *
                                    (r.air_mass_kg_m2 * r.tracer_mixing_ratio -
                                     l.air_mass_kg_m2 * l.tracer_mixing_ratio),
          .maximum_dissipation_speed_m_s = dissipation_speed,
          .maximum_wave_speed_m_s = maximum_wave_speed};
}
}  // namespace mps
