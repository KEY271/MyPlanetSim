#include "myplanetsim/dynamics/dry_hydrostatic_coupling.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
namespace mps {
DryHydrostaticCoupling couple_dry_hydrostatic_columns(
    const DryHydrostaticState& s, const DryHydrostaticDerived& d,
    const DryHydrostaticTransportTendency& h, const std::span<const Real> b,
    const Real g, const VerticalTransportScheme scheme,
    const VerticalLimiterKind limiter) {
  DryHydrostaticCoupling result;
  DryHydrostaticCouplingWorkspace workspace;
  couple_dry_hydrostatic_columns(s, d, h, b, g, scheme, limiter, result, workspace);
  return result;
}

void couple_dry_hydrostatic_columns(const DryHydrostaticState& s,
                                    const DryHydrostaticDerived& d,
                                    const DryHydrostaticTransportTendency& h,
                                    const std::span<const Real> b, const Real g,
                                    const VerticalTransportScheme scheme,
                                    const VerticalLimiterKind limiter,
                                    DryHydrostaticCoupling& out,
                                    DryHydrostaticCouplingWorkspace& workspace) {
  couple_dry_hydrostatic_columns(s, d, h, b, g, scheme, limiter, limiter, out,
                                 workspace);
}

void couple_dry_hydrostatic_columns(
    const DryHydrostaticState& s, const DryHydrostaticDerived& d,
    const DryHydrostaticTransportTendency& h, const std::span<const Real> b,
    const Real g, const VerticalTransportScheme scheme,
    const VerticalLimiterKind limiter,
    const VerticalLimiterKind potential_temperature_limiter,
    DryHydrostaticCoupling& out, DryHydrostaticCouplingWorkspace& workspace) {
  const auto n = d.cells * d.levels;
  if (h.air_mass.size() != n || h.momentum.size() != n ||
      h.potential_temperature_mass.size() != n || h.tracer_mass.size() != n ||
      b.size() != d.levels + 1)
    throw std::invalid_argument("dry coupling shape mismatch");
  out.surface_pressure_pa_s.resize(d.cells);
  out.interface_mass_flux_kg_m2_s.resize(d.cells * (d.levels + 1));
  out.tendency = h;
  out.maximum_continuity_residual_pa_s = 0.0;
  workspace.interface_coordinate.resize(d.levels + 1);
  workspace.center_coordinate.resize(d.levels);
  workspace.scalar.resize(d.levels);
  workspace.flux.resize(d.levels + 1);
  workspace.rhs.resize(d.levels);
  workspace.component.resize(d.levels);
  workspace.horizontal_component.resize(d.levels);
  for (std::size_t c = 0; c < d.cells; ++c) {
    const auto begin = c * d.levels;
    const std::span<const Real> horizontal_mass(h.air_mass.data() + begin, d.levels);
    const std::span<const Real> mass(d.air_mass_kg_m2.data() + begin, d.levels);
    diagnose_vertical_mass_flux(horizontal_mass, b, g, workspace.mass_flux);
    const auto& f = workspace.mass_flux;
    out.surface_pressure_pa_s[c] = f.surface_pressure_tendency_pa_s;
    out.maximum_continuity_residual_pa_s = std::max(
        out.maximum_continuity_residual_pa_s, std::abs(f.continuity_residual_pa_s));
    std::copy(f.interface_flux_kg_m2_s.begin(), f.interface_flux_kg_m2_s.end(),
              out.interface_mass_flux_kg_m2_s.begin() +
                  static_cast<std::ptrdiff_t>(c * (d.levels + 1)));
    workspace.interface_coordinate[0] = 0.0;
    for (std::size_t k = 0; k < d.levels; ++k) {
      workspace.interface_coordinate[k + 1] =
          workspace.interface_coordinate[k] + mass[k];
      workspace.center_coordinate[k] = 0.5 * (workspace.interface_coordinate[k] +
                                              workspace.interface_coordinate[k + 1]);
    }
    const auto apply = [&](const std::span<const Real> state,
                           const std::span<const Real> horizontal,
                           const VerticalLimiterKind field_limiter) {
      detail::vertical_scalar_rhs_unchecked(
          state, mass, horizontal, f, scheme, field_limiter,
          workspace.interface_coordinate, workspace.center_coordinate, workspace.scalar,
          workspace.flux, workspace.rhs);
    };
    apply(std::span<const Real>(s.potential_temperature_mass_k_kg_m2.data() + begin,
                                d.levels),
          std::span<const Real>(h.potential_temperature_mass.data() + begin, d.levels),
          potential_temperature_limiter);
    for (std::size_t k = 0; k < d.levels; ++k) {
      out.tendency.air_mass[begin + k] = f.target_air_mass_tendency_kg_m2_s[k];
      out.tendency.potential_temperature_mass[begin + k] = workspace.rhs[k];
    }
    apply(std::span<const Real>(s.tracer_mass_kg_m2.data() + begin, d.levels),
          std::span<const Real>(h.tracer_mass.data() + begin, d.levels), limiter);
    for (std::size_t k = 0; k < d.levels; ++k)
      out.tendency.tracer_mass[begin + k] = workspace.rhs[k];
    for (int component = 0; component < 3; ++component) {
      for (std::size_t k = 0; k < d.levels; ++k) {
        workspace.component[k] =
            s.horizontal_momentum_mass_kg_m_s[begin + k][component];
        workspace.horizontal_component[k] = h.momentum[begin + k][component];
      }
      apply(workspace.component, workspace.horizontal_component, limiter);
      for (std::size_t k = 0; k < d.levels; ++k)
        out.tendency.momentum[begin + k][component] = workspace.rhs[k];
    }
  }
}
}  // namespace mps
