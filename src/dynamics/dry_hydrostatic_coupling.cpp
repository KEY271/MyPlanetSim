#include "myplanetsim/dynamics/dry_hydrostatic_coupling.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#if defined(MPS_ENABLE_OPENMP)
#include <omp.h>
#endif

namespace mps {
namespace {

[[nodiscard]] std::size_t coupling_worker_count() {
#if defined(MPS_ENABLE_OPENMP)
  return static_cast<std::size_t>(std::max(1, omp_get_max_threads()));
#else
  return 1;
#endif
}

[[nodiscard]] std::size_t coupling_worker_index() {
#if defined(MPS_ENABLE_OPENMP)
  return static_cast<std::size_t>(omp_get_thread_num());
#else
  return 0;
#endif
}

}  // namespace

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
  if (d.tracer_count == 0 || s.tracer_count != d.tracer_count ||
      h.air_mass.size() != n || h.momentum.size() != n ||
      h.potential_temperature_mass.size() != n ||
      h.tracer_mass.size() != d.tracer_count * n || b.size() != d.levels + 1)
    throw std::invalid_argument("dry coupling shape mismatch");
  out.surface_pressure_pa_s.resize(d.cells);
  out.interface_mass_flux_kg_m2_s.resize(d.cells * (d.levels + 1));
  out.tendency = h;
  out.maximum_continuity_residual_pa_s = 0.0;
  workspace.columns.resize(coupling_worker_count());
  for (auto& column : workspace.columns) {
    column.interface_coordinate.resize(d.levels + 1);
    column.center_coordinate.resize(d.levels);
    column.scalar.resize(d.levels);
    column.flux.resize(d.levels + 1);
    column.rhs.resize(d.levels);
    column.component.resize(d.levels);
    column.horizontal_component.resize(d.levels);
  }
  workspace.continuity_residual_pa_s.resize(d.cells);
  workspace.failures.resize(d.cells);
#if defined(MPS_ENABLE_OPENMP)
#pragma omp parallel for schedule(static)
#endif
  for (std::ptrdiff_t cell_index = 0; cell_index < static_cast<std::ptrdiff_t>(d.cells);
       ++cell_index) {
    const auto c = static_cast<std::size_t>(cell_index);
    auto& column = workspace.columns[coupling_worker_index()];
    workspace.failures[c] = nullptr;
    try {
      const auto begin = c * d.levels;
      const std::span<const Real> horizontal_mass(h.air_mass.data() + begin, d.levels);
      const std::span<const Real> mass(d.air_mass_kg_m2.data() + begin, d.levels);
      diagnose_vertical_mass_flux(horizontal_mass, b, g, column.mass_flux);
      const auto& f = column.mass_flux;
      out.surface_pressure_pa_s[c] = f.surface_pressure_tendency_pa_s;
      workspace.continuity_residual_pa_s[c] = std::abs(f.continuity_residual_pa_s);
      std::copy(f.interface_flux_kg_m2_s.begin(), f.interface_flux_kg_m2_s.end(),
                out.interface_mass_flux_kg_m2_s.begin() +
                    static_cast<std::ptrdiff_t>(c * (d.levels + 1)));
      column.interface_coordinate[0] = 0.0;
      for (std::size_t k = 0; k < d.levels; ++k) {
        column.interface_coordinate[k + 1] = column.interface_coordinate[k] + mass[k];
        column.center_coordinate[k] =
            0.5 * (column.interface_coordinate[k] + column.interface_coordinate[k + 1]);
      }
      const auto apply = [&](const std::span<const Real> state,
                             const std::span<const Real> horizontal,
                             const VerticalLimiterKind field_limiter) {
        detail::vertical_scalar_rhs_unchecked(
            state, mass, horizontal, f, scheme, field_limiter,
            column.interface_coordinate, column.center_coordinate, column.scalar,
            column.flux, column.rhs);
      };
      apply(
          std::span<const Real>(s.potential_temperature_mass_k_kg_m2.data() + begin,
                                d.levels),
          std::span<const Real>(h.potential_temperature_mass.data() + begin, d.levels),
          potential_temperature_limiter);
      for (std::size_t k = 0; k < d.levels; ++k) {
        out.tendency.air_mass[begin + k] = f.target_air_mass_tendency_kg_m2_s[k];
        out.tendency.potential_temperature_mass[begin + k] = column.rhs[k];
      }
      for (std::size_t tracer = 0; tracer < d.tracer_count; ++tracer) {
        const auto tracer_begin = tracer * n + begin;
        const std::span<const Real> tracer_state(
            s.tracer_mass_kg_m2.data() + tracer_begin, d.levels);
        const std::span<const Real> horizontal_tracer(
            h.tracer_mass.data() + tracer_begin, d.levels);
        const bool tracer_is_inactive =
            std::ranges::all_of(tracer_state,
                                [](const Real value) { return value == 0.0; }) &&
            std::ranges::all_of(horizontal_tracer,
                                [](const Real value) { return value == 0.0; });
        if (tracer_is_inactive) {
          std::fill_n(out.tendency.tracer_mass.begin() +
                          static_cast<std::ptrdiff_t>(tracer_begin),
                      d.levels, 0.0);
        } else {
          apply(tracer_state, horizontal_tracer, limiter);
          for (std::size_t k = 0; k < d.levels; ++k)
            out.tendency.tracer_mass[tracer_begin + k] = column.rhs[k];
        }
      }
      for (int component = 0; component < 3; ++component) {
        for (std::size_t k = 0; k < d.levels; ++k) {
          column.component[k] = s.horizontal_momentum_mass_kg_m_s[begin + k][component];
          column.horizontal_component[k] = h.momentum[begin + k][component];
        }
        apply(column.component, column.horizontal_component, limiter);
        for (std::size_t k = 0; k < d.levels; ++k)
          out.tendency.momentum[begin + k][component] = column.rhs[k];
      }
    } catch (...) {
      workspace.failures[c] = std::current_exception();
    }
  }
  for (std::size_t c = 0; c < d.cells; ++c) {
    if (workspace.failures[c] != nullptr) std::rethrow_exception(workspace.failures[c]);
    out.maximum_continuity_residual_pa_s = std::max(
        out.maximum_continuity_residual_pa_s, workspace.continuity_residual_pa_s[c]);
  }
}
}  // namespace mps
