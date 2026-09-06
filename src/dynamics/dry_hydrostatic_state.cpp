#include "myplanetsim/dynamics/dry_hydrostatic_state.hpp"

#include <cmath>
#include <stdexcept>

namespace mps {
namespace {
// Second-order face reconstruction lets a mixing ratio that is physically zero pick up
// a vanishingly small negative value from flux cancellation. Rejecting those is a false
// alarm, so nonnegativity is checked to roundoff, the way tangency already is. A tracer
// that is genuinely being driven negative fails this by many orders of magnitude.
constexpr Real kTracerRoundoffTolerance = 1e-12;

void require_shape(const DryHydrostaticState& state, const std::size_t levels) {
  const auto cells = state.surface_pressure_pa.size();
  const auto volume = cells * levels;
  if (cells == 0 || levels == 0 ||
      state.horizontal_momentum_mass_kg_m_s.size() != volume ||
      state.potential_temperature_mass_k_kg_m2.size() != volume ||
      state.tracer_mass_kg_m2.size() != volume ||
      (!state.surface_temperature_k.empty() &&
       state.surface_temperature_k.size() != cells)) {
    throw std::invalid_argument("dry hydrostatic state shape is invalid");
  }
}
}  // namespace

std::vector<Real> flatten_dry_hydrostatic_state(const DryHydrostaticState& state,
                                                const std::size_t levels) {
  require_shape(state, levels);
  const auto cells = state.surface_pressure_pa.size();
  const auto volume = cells * levels;
  std::vector<Real> values(cells + 5 * volume);
  std::copy(state.surface_pressure_pa.begin(), state.surface_pressure_pa.end(),
            values.begin());
  std::size_t cursor = cells;
  for (const auto v : state.horizontal_momentum_mass_kg_m_s) {
    values[cursor++] = v.x;
    values[cursor++] = v.y;
    values[cursor++] = v.z;
  }
  std::copy(state.potential_temperature_mass_k_kg_m2.begin(),
            state.potential_temperature_mass_k_kg_m2.end(),
            values.begin() + static_cast<std::ptrdiff_t>(cursor));
  cursor += volume;
  std::copy(state.tracer_mass_kg_m2.begin(), state.tracer_mass_kg_m2.end(),
            values.begin() + static_cast<std::ptrdiff_t>(cursor));
  return values;
}

DryHydrostaticState unflatten_dry_hydrostatic_state(const Real time_s,
                                                    const std::uint64_t step,
                                                    const std::span<const Real> values,
                                                    const std::size_t cells,
                                                    const std::size_t levels) {
  const auto volume = cells * levels;
  if (cells == 0 || levels == 0 || values.size() != cells + 5 * volume)
    throw std::invalid_argument("dry hydrostatic flat state size does not match shape");
  DryHydrostaticState state{.time_s = time_s,
                            .step = step,
                            .surface_pressure_pa = {},
                            .horizontal_momentum_mass_kg_m_s = {},
                            .potential_temperature_mass_k_kg_m2 = {},
                            .tracer_mass_kg_m2 = {},
                            .surface_temperature_k = {}};
  state.surface_pressure_pa.assign(values.begin(),
                                   values.begin() + static_cast<std::ptrdiff_t>(cells));
  state.horizontal_momentum_mass_kg_m_s.resize(volume);
  std::size_t cursor = cells;
  for (auto& v : state.horizontal_momentum_mass_kg_m_s)
    v = {values[cursor++], values[cursor++], values[cursor++]};
  state.potential_temperature_mass_k_kg_m2.assign(
      values.begin() + static_cast<std::ptrdiff_t>(cursor),
      values.begin() + static_cast<std::ptrdiff_t>(cursor + volume));
  cursor += volume;
  state.tracer_mass_kg_m2.assign(values.begin() + static_cast<std::ptrdiff_t>(cursor),
                                 values.end());
  return state;
}

std::vector<Real> flatten_dry_hydrostatic_surface_state(
    const DryHydrostaticState& state, const std::size_t levels) {
  if (state.surface_temperature_k.size() != state.surface_pressure_pa.size())
    throw std::invalid_argument("surface state must contain one temperature per cell");
  auto values = flatten_dry_hydrostatic_state(state, levels);
  values.insert(values.end(), state.surface_temperature_k.begin(),
                state.surface_temperature_k.end());
  return values;
}

DryHydrostaticState unflatten_dry_hydrostatic_surface_state(
    const Real time_s, const std::uint64_t step, const std::span<const Real> values,
    const std::size_t cells, const std::size_t levels) {
  const std::size_t atmospheric_size = cells + 5 * cells * levels;
  if (cells == 0 || values.size() != atmospheric_size + cells)
    throw std::invalid_argument(
        "dry hydrostatic surface flat state size does not match shape");
  auto state = unflatten_dry_hydrostatic_state(
      time_s, step, values.first(atmospheric_size), cells, levels);
  state.surface_temperature_k.assign(
      values.begin() + static_cast<std::ptrdiff_t>(atmospheric_size), values.end());
  return state;
}

DryHydrostaticDerived diagnose_dry_hydrostatic_state(
    const DryHydrostaticState& state, const AtmosphericHybridCoordinate& coordinate,
    const PlanetParameters& planet) {
  return diagnose_dry_hydrostatic_state(state, coordinate, planet, {});
}

DryHydrostaticDerived diagnose_dry_hydrostatic_state(
    const DryHydrostaticState& state, const AtmosphericHybridCoordinate& coordinate,
    const PlanetParameters& planet,
    const std::span<const Real> surface_geopotential_m2_s2) {
  DryHydrostaticDerived out;
  HybridPressureGeometry geometry;
  HydrostaticColumn column;
  std::vector<Real> theta(coordinate.levels());
  diagnose_dry_hydrostatic_state(state, coordinate, planet, surface_geopotential_m2_s2,
                                 out, geometry, column, theta);
  return out;
}

void diagnose_dry_hydrostatic_state(
    const DryHydrostaticState& state, const AtmosphericHybridCoordinate& coordinate,
    const PlanetParameters& planet,
    const std::span<const Real> surface_geopotential_m2_s2, DryHydrostaticDerived& out,
    HybridPressureGeometry& geometry, HydrostaticColumn& hydro,
    const std::span<Real> theta) {
  require_shape(state, coordinate.levels());
  if (!surface_geopotential_m2_s2.empty() &&
      surface_geopotential_m2_s2.size() != state.surface_pressure_pa.size())
    throw std::invalid_argument("surface orography shape does not match state");
  out.cells = state.surface_pressure_pa.size();
  out.levels = coordinate.levels();
  if (theta.size() != out.levels)
    throw std::invalid_argument("dry diagnostic theta workspace shape mismatch");
  const auto volume = out.cells * out.levels;
  out.pressure_pa.resize(volume);
  out.exner_half.resize(out.cells * (out.levels + 1));
  out.exner_full.resize(volume);
  out.air_mass_kg_m2.resize(volume);
  out.velocity_m_s.resize(volume);
  out.potential_temperature_k.resize(volume);
  out.tracer_mixing_ratio.resize(volume);
  out.temperature_k.resize(volume);
  out.geopotential_m2_s2.resize(volume);
  for (std::size_t c = 0; c < out.cells; ++c) {
    coordinate.geometry(state.surface_pressure_pa[c], planet.gravity_m_s2,
                        planet.gas_constant_j_kg_k, planet.heat_capacity_cp_j_kg_k,
                        planet.reference_pressure_pa, geometry);
    std::copy(
        geometry.exner_half.begin(), geometry.exner_half.end(),
        out.exner_half.begin() + static_cast<std::ptrdiff_t>(c * (out.levels + 1)));
    for (std::size_t k = 0; k < out.levels; ++k) {
      const auto n = dry_hydrostatic_offset(c, k, out.levels);
      const auto mass = geometry.air_mass_kg_m2[k];
      out.pressure_pa[n] = geometry.pressure_full_pa[k];
      out.exner_full[n] = geometry.exner_full[k];
      out.air_mass_kg_m2[n] = mass;
      out.velocity_m_s[n] = state.horizontal_momentum_mass_kg_m_s[n] / mass;
      theta[k] = state.potential_temperature_mass_k_kg_m2[n] / mass;
      out.potential_temperature_k[n] = theta[k];
      out.tracer_mixing_ratio[n] = state.tracer_mass_kg_m2[n] / mass;
      out.temperature_k[n] = theta[k] * geometry.exner_full[k];
    }
    const Real surface_geopotential =
        surface_geopotential_m2_s2.empty() ? 0.0 : surface_geopotential_m2_s2[c];
    integrate_hydrostatic_column(geometry, theta, planet.heat_capacity_cp_j_kg_k,
                                 planet.gravity_m_s2, surface_geopotential, hydro);
    for (std::size_t k = 0; k < out.levels; ++k)
      out.geopotential_m2_s2[dry_hydrostatic_offset(c, k, out.levels)] =
          hydro.geopotential_full_m2_s2[k];
  }
}

void validate_dry_hydrostatic_state(const DryHydrostaticState& state,
                                    const DryHydrostaticDerived& d,
                                    const std::span<const Vec3> centres,
                                    const Real min_ps, const Real max_ps,
                                    const Real floor, const bool nonnegative) {
  require_shape(state, d.levels);
  if (d.cells != state.surface_pressure_pa.size() || centres.size() != d.cells ||
      d.air_mass_kg_m2.size() != d.cells * d.levels)
    throw std::invalid_argument("dry hydrostatic derived shape is invalid");
  if (!std::isfinite(state.time_s) || state.time_s < 0.0)
    throw std::runtime_error("dry hydrostatic time is invalid");
  for (std::size_t c = 0; c < d.cells; ++c) {
    if (!std::isfinite(state.surface_pressure_pa[c]) ||
        state.surface_pressure_pa[c] < min_ps || state.surface_pressure_pa[c] > max_ps)
      throw std::runtime_error("dry hydrostatic surface pressure is invalid");
    for (std::size_t k = 0; k < d.levels; ++k) {
      const auto n = dry_hydrostatic_offset(c, k, d.levels);
      if (!(d.air_mass_kg_m2[n] > 0.0) || !std::isfinite(d.temperature_k[n]) ||
          d.temperature_k[n] < floor ||
          !is_finite(state.horizontal_momentum_mass_kg_m_s[n]) ||
          std::abs(dot(state.horizontal_momentum_mass_kg_m_s[n], centres[c])) >
              1e-10 * std::max(1.0, norm(state.horizontal_momentum_mass_kg_m_s[n])) ||
          !std::isfinite(d.potential_temperature_k[n]) ||
          !(d.potential_temperature_k[n] > 0.0) ||
          !std::isfinite(d.tracer_mixing_ratio[n]) ||
          (nonnegative && d.tracer_mixing_ratio[n] < -kTracerRoundoffTolerance))
        throw std::runtime_error("dry hydrostatic cell-layer invariant failed");
    }
  }
  for (const Real temperature : state.surface_temperature_k)
    if (!(temperature > 0.0) || !std::isfinite(temperature))
      throw std::runtime_error("dry hydrostatic surface temperature is invalid");
}
}  // namespace mps
