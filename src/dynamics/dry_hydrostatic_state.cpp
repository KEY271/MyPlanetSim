#include "myplanetsim/dynamics/dry_hydrostatic_state.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

#if defined(MPS_ENABLE_OPENMP)
#include <omp.h>
#endif

namespace mps {
namespace {
// Second-order face reconstruction lets a mixing ratio that is physically zero pick up
// a vanishingly small negative value from flux cancellation. Rejecting those is a false
// alarm, so nonnegativity is checked to roundoff, the way tangency already is. A tracer
// that is genuinely being driven negative fails this by many orders of magnitude.
constexpr Real kTracerRoundoffTolerance = 1e-12;

[[nodiscard]] std::size_t diagnosis_worker_count() {
#if defined(MPS_ENABLE_OPENMP)
  return static_cast<std::size_t>(std::max(1, omp_get_max_threads()));
#else
  return 1;
#endif
}

[[nodiscard]] std::size_t diagnosis_worker_index() {
#if defined(MPS_ENABLE_OPENMP)
  return static_cast<std::size_t>(omp_get_thread_num());
#else
  return 0;
#endif
}

void require_shape(const DryHydrostaticState& state, const std::size_t levels) {
  const auto cells = state.surface_pressure_pa.size();
  const auto volume = cells * levels;
  if (cells == 0 || levels == 0 ||
      state.horizontal_momentum_mass_kg_m_s.size() != volume ||
      state.potential_temperature_mass_k_kg_m2.size() != volume ||
      state.tracer_count == 0 ||
      state.tracer_mass_kg_m2.size() != state.tracer_count * volume ||
      (!state.surface_temperature_k.empty() &&
       state.surface_temperature_k.size() != cells)) {
    throw std::invalid_argument("dry hydrostatic state shape is invalid");
  }
}

void require_valid_moist_ledgers(const DryHydrostaticState& state) {
  const Real ledgers[]{
      state.cumulative_convective_precipitation_kg,
      state.cumulative_grid_scale_precipitation_kg,
      state.cumulative_evaporation_kg,
      state.cumulative_runoff_kg,
      state.cumulative_ocean_water_change_kg,
      state.cumulative_external_outflow_kg,
  };
  for (const Real value : ledgers)
    if (!std::isfinite(value))
      throw std::invalid_argument("moist checkpoint contains a non-finite ledger");
}
}  // namespace

std::span<Real> dry_hydrostatic_tracer_component(DryHydrostaticState& state,
                                                 const std::size_t tracer,
                                                 const std::size_t cells,
                                                 const std::size_t levels) {
  const auto volume = cells * levels;
  if (tracer >= state.tracer_count ||
      state.tracer_mass_kg_m2.size() != state.tracer_count * volume)
    throw std::out_of_range("dry hydrostatic tracer component is out of range");
  return std::span<Real>(state.tracer_mass_kg_m2).subspan(tracer * volume, volume);
}

std::span<const Real> dry_hydrostatic_tracer_component(const DryHydrostaticState& state,
                                                       const std::size_t tracer,
                                                       const std::size_t cells,
                                                       const std::size_t levels) {
  const auto volume = cells * levels;
  if (tracer >= state.tracer_count ||
      state.tracer_mass_kg_m2.size() != state.tracer_count * volume)
    throw std::out_of_range("dry hydrostatic tracer component is out of range");
  return std::span<const Real>(state.tracer_mass_kg_m2)
      .subspan(tracer * volume, volume);
}

std::vector<Real> flatten_dry_hydrostatic_state(const DryHydrostaticState& state,
                                                const std::size_t levels) {
  require_shape(state, levels);
  if (state.tracer_count != 1)
    throw std::invalid_argument("legacy dry checkpoint requires one tracer");
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
                            .tracer_count = 1,
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

std::string dry_hydrostatic_multitracer_checkpoint_layout(
    const TracerRegistry& registry, const bool include_surface) {
  std::string result(include_surface ? kDryHydrostaticMultitracerSurfaceCheckpointLayout
                                     : kDryHydrostaticMultitracerCheckpointLayout);
  result += '|';
  for (std::size_t index = 0; index < registry.size(); ++index) {
    if (index != 0) result += ',';
    result += registry[index].name;
    result += ':';
    result += tracer_role_name(registry[index].role);
  }
  return result;
}

std::vector<Real> flatten_dry_hydrostatic_multitracer_state(
    const DryHydrostaticState& state, const std::size_t levels,
    const bool include_surface) {
  require_shape(state, levels);
  const auto cells = state.surface_pressure_pa.size();
  if (include_surface && state.surface_temperature_k.size() != cells)
    throw std::invalid_argument(
        "multi-tracer surface checkpoint requires one temperature per cell");
  if (!include_surface && !state.surface_temperature_k.empty())
    throw std::invalid_argument(
        "atmosphere-only multi-tracer checkpoint rejects a surface state");
  const auto volume = cells * levels;
  std::vector<Real> values(cells + (4 + state.tracer_count) * volume +
                           (include_surface ? cells : 0));
  std::copy(state.surface_pressure_pa.begin(), state.surface_pressure_pa.end(),
            values.begin());
  std::size_t cursor = cells;
  for (const auto value : state.horizontal_momentum_mass_kg_m_s) {
    values[cursor++] = value.x;
    values[cursor++] = value.y;
    values[cursor++] = value.z;
  }
  std::copy(state.potential_temperature_mass_k_kg_m2.begin(),
            state.potential_temperature_mass_k_kg_m2.end(),
            values.begin() + static_cast<std::ptrdiff_t>(cursor));
  cursor += volume;
  std::copy(state.tracer_mass_kg_m2.begin(), state.tracer_mass_kg_m2.end(),
            values.begin() + static_cast<std::ptrdiff_t>(cursor));
  cursor += state.tracer_count * volume;
  if (include_surface)
    std::copy(state.surface_temperature_k.begin(), state.surface_temperature_k.end(),
              values.begin() + static_cast<std::ptrdiff_t>(cursor));
  return values;
}

DryHydrostaticState unflatten_dry_hydrostatic_multitracer_state(
    const Real time_s, const std::uint64_t step, const std::span<const Real> values,
    const std::size_t cells, const std::size_t levels, const std::size_t tracer_count,
    const bool include_surface) {
  const auto volume = cells * levels;
  const auto expected =
      cells + (4 + tracer_count) * volume + (include_surface ? cells : 0);
  if (cells == 0 || levels == 0 || tracer_count == 0 || values.size() != expected)
    throw std::invalid_argument(
        "multi-tracer flat state size does not match registry and shape");
  DryHydrostaticState state{
      .time_s = time_s, .step = step, .tracer_count = tracer_count};
  state.surface_pressure_pa.assign(values.begin(),
                                   values.begin() + static_cast<std::ptrdiff_t>(cells));
  state.horizontal_momentum_mass_kg_m_s.resize(volume);
  std::size_t cursor = cells;
  for (auto& value : state.horizontal_momentum_mass_kg_m_s)
    value = {values[cursor++], values[cursor++], values[cursor++]};
  state.potential_temperature_mass_k_kg_m2.assign(
      values.begin() + static_cast<std::ptrdiff_t>(cursor),
      values.begin() + static_cast<std::ptrdiff_t>(cursor + volume));
  cursor += volume;
  state.tracer_mass_kg_m2.assign(
      values.begin() + static_cast<std::ptrdiff_t>(cursor),
      values.begin() + static_cast<std::ptrdiff_t>(cursor + tracer_count * volume));
  cursor += tracer_count * volume;
  if (include_surface)
    state.surface_temperature_k.assign(
        values.begin() + static_cast<std::ptrdiff_t>(cursor), values.end());
  require_shape(state, levels);
  return state;
}

std::string dry_hydrostatic_moist_checkpoint_layout(const TracerRegistry& registry) {
  std::string result(kDryHydrostaticMoistCheckpointLayout);
  result +=
      dry_hydrostatic_multitracer_checkpoint_layout(registry, true)
          .substr(
              std::string(kDryHydrostaticMultitracerSurfaceCheckpointLayout).size());
  return result;
}

std::vector<Real> flatten_dry_hydrostatic_moist_state(const DryHydrostaticState& state,
                                                      const std::size_t levels) {
  const auto cells = state.surface_pressure_pa.size();
  if (state.land_water_kg_m2.size() != cells)
    throw std::invalid_argument("moist checkpoint requires one land bucket per cell");
  for (const Real water : state.land_water_kg_m2)
    if (water < 0.0 || !std::isfinite(water))
      throw std::invalid_argument("moist checkpoint contains invalid land water");
  require_valid_moist_ledgers(state);
  auto values = flatten_dry_hydrostatic_multitracer_state(state, levels, true);
  values.insert(values.end(), state.land_water_kg_m2.begin(),
                state.land_water_kg_m2.end());
  values.push_back(state.cumulative_convective_precipitation_kg);
  values.push_back(state.cumulative_grid_scale_precipitation_kg);
  values.push_back(state.cumulative_evaporation_kg);
  values.push_back(state.cumulative_runoff_kg);
  values.push_back(state.cumulative_ocean_water_change_kg);
  values.push_back(state.cumulative_external_outflow_kg);
  return values;
}

DryHydrostaticState unflatten_dry_hydrostatic_moist_state(
    const Real time_s, const std::uint64_t step, const std::span<const Real> values,
    const std::size_t cells, const std::size_t levels, const std::size_t tracer_count) {
  const auto atmospheric_size = cells + (4 + tracer_count) * cells * levels + cells;
  if (values.size() != atmospheric_size + cells + 6)
    throw std::invalid_argument("moist checkpoint size does not match state shape");
  auto state = unflatten_dry_hydrostatic_multitracer_state(
      time_s, step, values.first(atmospheric_size), cells, levels, tracer_count, true);
  std::size_t cursor = atmospheric_size;
  state.land_water_kg_m2.assign(
      values.begin() + static_cast<std::ptrdiff_t>(cursor),
      values.begin() + static_cast<std::ptrdiff_t>(cursor + cells));
  cursor += cells;
  state.cumulative_convective_precipitation_kg = values[cursor++];
  state.cumulative_grid_scale_precipitation_kg = values[cursor++];
  state.cumulative_evaporation_kg = values[cursor++];
  state.cumulative_runoff_kg = values[cursor++];
  state.cumulative_ocean_water_change_kg = values[cursor++];
  state.cumulative_external_outflow_kg = values[cursor];
  for (const Real water : state.land_water_kg_m2)
    if (water < 0.0 || !std::isfinite(water))
      throw std::invalid_argument("moist checkpoint contains invalid land water");
  require_valid_moist_ledgers(state);
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
  out.tracer_count = state.tracer_count;
  if (theta.size() != out.levels)
    throw std::invalid_argument("dry diagnostic theta workspace shape mismatch");
  const auto volume = out.cells * out.levels;
  out.pressure_pa.resize(volume);
  out.exner_half.resize(out.cells * (out.levels + 1));
  out.exner_full.resize(volume);
  out.air_mass_kg_m2.resize(volume);
  out.velocity_m_s.resize(volume);
  out.potential_temperature_k.resize(volume);
  out.tracer_mixing_ratio.resize(state.tracer_count * volume);
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
      for (std::size_t tracer = 0; tracer < state.tracer_count; ++tracer) {
        const auto q =
            dry_hydrostatic_tracer_offset(tracer, c, k, out.cells, out.levels);
        out.tracer_mixing_ratio[q] = state.tracer_mass_kg_m2[q] / mass;
      }
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

void diagnose_dry_hydrostatic_state(
    const DryHydrostaticState& state, const AtmosphericHybridCoordinate& coordinate,
    const PlanetParameters& planet,
    const std::span<const Real> surface_geopotential_m2_s2, DryHydrostaticDerived& out,
    DryHydrostaticDiagnosisWorkspace& workspace) {
  require_shape(state, coordinate.levels());
  if (!surface_geopotential_m2_s2.empty() &&
      surface_geopotential_m2_s2.size() != state.surface_pressure_pa.size())
    throw std::invalid_argument("surface orography shape does not match state");
  out.cells = state.surface_pressure_pa.size();
  out.levels = coordinate.levels();
  out.tracer_count = state.tracer_count;
  const auto volume = out.cells * out.levels;
  out.pressure_pa.resize(volume);
  out.exner_half.resize(out.cells * (out.levels + 1));
  out.exner_full.resize(volume);
  out.air_mass_kg_m2.resize(volume);
  out.velocity_m_s.resize(volume);
  out.potential_temperature_k.resize(volume);
  out.tracer_mixing_ratio.resize(state.tracer_count * volume);
  out.temperature_k.resize(volume);
  out.geopotential_m2_s2.resize(volume);
  workspace.columns.resize(diagnosis_worker_count());
  for (auto& column : workspace.columns)
    column.potential_temperature.resize(out.levels);
  workspace.failures.resize(out.cells);
#if defined(MPS_ENABLE_OPENMP)
#pragma omp parallel for schedule(static)
#endif
  for (std::ptrdiff_t cell_index = 0;
       cell_index < static_cast<std::ptrdiff_t>(out.cells); ++cell_index) {
    const auto cell = static_cast<std::size_t>(cell_index);
    auto& column = workspace.columns[diagnosis_worker_index()];
    workspace.failures[cell] = nullptr;
    try {
      coordinate.geometry(state.surface_pressure_pa[cell], planet.gravity_m_s2,
                          planet.gas_constant_j_kg_k, planet.heat_capacity_cp_j_kg_k,
                          planet.reference_pressure_pa, column.geometry);
      std::copy(column.geometry.exner_half.begin(), column.geometry.exner_half.end(),
                out.exner_half.begin() +
                    static_cast<std::ptrdiff_t>(cell * (out.levels + 1)));
      for (std::size_t level = 0; level < out.levels; ++level) {
        const auto n = dry_hydrostatic_offset(cell, level, out.levels);
        const auto mass = column.geometry.air_mass_kg_m2[level];
        out.pressure_pa[n] = column.geometry.pressure_full_pa[level];
        out.exner_full[n] = column.geometry.exner_full[level];
        out.air_mass_kg_m2[n] = mass;
        out.velocity_m_s[n] = state.horizontal_momentum_mass_kg_m_s[n] / mass;
        column.potential_temperature[level] =
            state.potential_temperature_mass_k_kg_m2[n] / mass;
        out.potential_temperature_k[n] = column.potential_temperature[level];
        for (std::size_t tracer = 0; tracer < state.tracer_count; ++tracer) {
          const auto q =
              dry_hydrostatic_tracer_offset(tracer, cell, level, out.cells, out.levels);
          out.tracer_mixing_ratio[q] = state.tracer_mass_kg_m2[q] / mass;
        }
        out.temperature_k[n] =
            column.potential_temperature[level] * column.geometry.exner_full[level];
      }
      const Real surface_geopotential =
          surface_geopotential_m2_s2.empty() ? 0.0 : surface_geopotential_m2_s2[cell];
      integrate_hydrostatic_column(column.geometry, column.potential_temperature,
                                   planet.heat_capacity_cp_j_kg_k, planet.gravity_m_s2,
                                   surface_geopotential, column.hydrostatic);
      for (std::size_t level = 0; level < out.levels; ++level)
        out.geopotential_m2_s2[dry_hydrostatic_offset(cell, level, out.levels)] =
            column.hydrostatic.geopotential_full_m2_s2[level];
    } catch (...) {
      workspace.failures[cell] = std::current_exception();
    }
  }
  for (std::size_t cell = 0; cell < out.cells; ++cell)
    if (workspace.failures[cell] != nullptr)
      std::rethrow_exception(workspace.failures[cell]);
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
          !(d.potential_temperature_k[n] > 0.0))
        throw std::runtime_error("dry hydrostatic cell-layer invariant failed");
      for (std::size_t tracer = 0; tracer < state.tracer_count; ++tracer) {
        const auto q = dry_hydrostatic_tracer_offset(tracer, c, k, d.cells, d.levels);
        if (!std::isfinite(d.tracer_mixing_ratio[q]) ||
            (nonnegative && d.tracer_mixing_ratio[q] < -kTracerRoundoffTolerance))
          throw std::runtime_error("dry hydrostatic tracer invariant failed");
      }
    }
  }
  for (const Real temperature : state.surface_temperature_k)
    if (!(temperature > 0.0) || !std::isfinite(temperature))
      throw std::runtime_error("dry hydrostatic surface temperature is invalid");
  if (!state.land_water_kg_m2.empty() && state.land_water_kg_m2.size() != d.cells)
    throw std::runtime_error("dry hydrostatic land-water shape is invalid");
  for (const Real water : state.land_water_kg_m2)
    if (water < 0.0 || !std::isfinite(water))
      throw std::runtime_error("dry hydrostatic land water is invalid");
  for (const Real ledger :
       {state.cumulative_convective_precipitation_kg,
        state.cumulative_grid_scale_precipitation_kg, state.cumulative_evaporation_kg,
        state.cumulative_runoff_kg, state.cumulative_ocean_water_change_kg,
        state.cumulative_external_outflow_kg})
    if (!std::isfinite(ledger))
      throw std::runtime_error("dry hydrostatic water ledger is invalid");
}
}  // namespace mps
