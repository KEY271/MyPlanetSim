#include "myplanetsim/dynamics/vertical_column_driver.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <numeric>
#include <stdexcept>

#include "myplanetsim/thermodynamics/dry_thermodynamics.hpp"

namespace mps {
namespace {
struct VerticalRhs {
  Real surface_pressure_pa_s;
  std::vector<Real> horizontal_air_mass_kg_m2_s;
  std::vector<Real> horizontal_theta_mass_k_kg_m2_s;
  std::vector<Real> horizontal_tracer_mass_kg_m2_s;
  std::vector<Real> theta_mass_k_kg_m2_s;
  std::vector<Real> tracer_mass_kg_m2_s;
  VerticalMassFlux mass_flux;
};

VerticalRhs rhs(const ExperimentConfig& config,
                const AtmosphericHybridCoordinate& coordinate,
                const VerticalColumnState& state) {
  const auto geometry = coordinate.geometry(
      state.surface_pressure_pa, config.planet.gravity_m_s2,
      config.planet.gas_constant_j_kg_k, config.planet.heat_capacity_cp_j_kg_k,
      config.planet.reference_pressure_pa);
  const std::size_t nz = coordinate.levels();
  std::vector<Real> h(nz, 0.0);
  if (config.vertical.test_case == VerticalTestCase::kMovingSurfacePressure ||
      config.vertical.test_case == VerticalTestCase::kManufacturedTransport) {
    const Real phase = 2.0 * std::numbers::pi *
                       (state.time_s - config.run.start_time_s) /
                       (config.run.end_time_s - config.run.start_time_s);
    const Real forcing = config.vertical.forcing_amplitude * std::cos(phase) /
                         config.planet.gravity_m_s2;
    for (std::size_t k = 0; k < nz; ++k) {
      const Real centered_index =
          static_cast<Real>(k) - 0.5 * static_cast<Real>(nz - 1);
      const Real weight = 1.0 + 0.25 * centered_index / static_cast<Real>(nz);
      h[k] = forcing * weight / static_cast<Real>(nz);
    }
  }
  const auto flux = diagnose_vertical_mass_flux(h, coordinate.coefficients().b_half,
                                                config.planet.gravity_m_s2);
  std::vector<Real> h_theta(nz);
  std::vector<Real> h_tracer(nz);
  for (std::size_t k = 0; k < nz; ++k) {
    const Real theta =
        state.potential_temperature_mass_k_kg_m2[k] / geometry.air_mass_kg_m2[k];
    const Real tracer = state.tracer_mass_kg_m2[k] / geometry.air_mass_kg_m2[k];
    h_theta[k] = theta * h[k];
    h_tracer[k] = tracer * h[k];
  }
  return {.surface_pressure_pa_s = flux.surface_pressure_tendency_pa_s,
          .horizontal_air_mass_kg_m2_s = h,
          .horizontal_theta_mass_k_kg_m2_s = h_theta,
          .horizontal_tracer_mass_kg_m2_s = h_tracer,
          .theta_mass_k_kg_m2_s = vertical_scalar_rhs(
              state.potential_temperature_mass_k_kg_m2, geometry.air_mass_kg_m2,
              h_theta, flux, config.vertical.transport_scheme, config.vertical.limiter),
          .tracer_mass_kg_m2_s = vertical_scalar_rhs(
              state.tracer_mass_kg_m2, geometry.air_mass_kg_m2, h_tracer, flux,
              config.vertical.transport_scheme, config.vertical.limiter),
          .mass_flux = flux};
}

VerticalColumnState plus(const VerticalColumnState& state, const VerticalRhs& tendency,
                         const Real factor) {
  VerticalColumnState result = state;
  result.surface_pressure_pa += factor * tendency.surface_pressure_pa_s;
  for (std::size_t k = 0; k < result.potential_temperature_mass_k_kg_m2.size(); ++k) {
    result.potential_temperature_mass_k_kg_m2[k] +=
        factor * tendency.theta_mass_k_kg_m2_s[k];
    result.tracer_mass_kg_m2[k] += factor * tendency.tracer_mass_kg_m2_s[k];
  }
  return result;
}

VerticalColumnState blend(const VerticalColumnState& first, const Real first_weight,
                          const VerticalColumnState& second, const Real second_weight) {
  VerticalColumnState result = first;
  result.surface_pressure_pa = first_weight * first.surface_pressure_pa +
                               second_weight * second.surface_pressure_pa;
  for (std::size_t k = 0; k < result.potential_temperature_mass_k_kg_m2.size(); ++k) {
    result.potential_temperature_mass_k_kg_m2[k] =
        first_weight * first.potential_temperature_mass_k_kg_m2[k] +
        second_weight * second.potential_temperature_mass_k_kg_m2[k];
    result.tracer_mass_kg_m2[k] = first_weight * first.tracer_mass_kg_m2[k] +
                                  second_weight * second.tracer_mass_kg_m2[k];
  }
  return result;
}

void validate(const ExperimentConfig& config,
              const AtmosphericHybridCoordinate& coordinate,
              const VerticalColumnState& state) {
  if (state.time_s < config.run.start_time_s || state.time_s > config.run.end_time_s) {
    throw std::invalid_argument("vertical column time is outside the configured run");
  }
  const auto geometry = coordinate.geometry(
      state.surface_pressure_pa, config.planet.gravity_m_s2,
      config.planet.gas_constant_j_kg_k, config.planet.heat_capacity_cp_j_kg_k,
      config.planet.reference_pressure_pa);
  validate_vertical_column_state(state, geometry.air_mass_kg_m2,
                                 config.vertical.temperature_floor_k,
                                 geometry.exner_full);
}

[[nodiscard]] Real sum(const std::vector<Real>& values) {
  return std::accumulate(values.begin(), values.end(), 0.0);
}

[[nodiscard]] bool pressure_is_in_range(const ExperimentConfig& config,
                                        const Real surface_pressure_pa) {
  return std::isfinite(surface_pressure_pa) &&
         surface_pressure_pa >= config.vertical.minimum_surface_pressure_pa &&
         surface_pressure_pa <= config.vertical.maximum_surface_pressure_pa;
}

[[nodiscard]] Real pressure_limited_time_step(const ExperimentConfig& config,
                                              const Real surface_pressure_pa,
                                              const Real surface_pressure_tendency_pa_s,
                                              const Real maximum_time_step_s) {
  Real available_pressure_pa = std::numeric_limits<Real>::infinity();
  if (surface_pressure_tendency_pa_s > 0.0) {
    available_pressure_pa =
        config.vertical.maximum_surface_pressure_pa - surface_pressure_pa;
  } else if (surface_pressure_tendency_pa_s < 0.0) {
    available_pressure_pa =
        surface_pressure_pa - config.vertical.minimum_surface_pressure_pa;
  } else {
    return maximum_time_step_s;
  }
  if (!(available_pressure_pa > 0.0)) {
    throw std::runtime_error(
        "vertical forcing drives surface pressure outside the configured range");
  }
  const Real pressure_limit =
      available_pressure_pa / std::abs(surface_pressure_tendency_pa_s);
  if (pressure_limit < maximum_time_step_s) {
    return std::nextafter(pressure_limit, 0.0);
  }
  return maximum_time_step_s;
}

void halve_time_step(const VerticalColumnState& state, Real& time_step_s) {
  time_step_s *= 0.5;
  if (!(time_step_s > 0.0) || state.time_s + time_step_s == state.time_s) {
    throw std::runtime_error(
        "surface-pressure or stage-CFL constraint produced a zero time step");
  }
}

void accumulate_budget(VerticalColumnBudget& budget, const VerticalRhs& first,
                       const VerticalRhs& second, const VerticalRhs& third,
                       const Real time_step_s) {
  const auto rk_integral = [time_step_s](const std::vector<Real>& first_values,
                                         const std::vector<Real>& second_values,
                                         const std::vector<Real>& third_values) {
    return time_step_s * (sum(first_values) / 6.0 + sum(second_values) / 6.0 +
                          2.0 * sum(third_values) / 3.0);
  };
  budget.integrated_dry_mass_source_kg_m2 +=
      rk_integral(first.horizontal_air_mass_kg_m2_s, second.horizontal_air_mass_kg_m2_s,
                  third.horizontal_air_mass_kg_m2_s);
  budget.integrated_potential_temperature_mass_source_k_kg_m2 += rk_integral(
      first.horizontal_theta_mass_k_kg_m2_s, second.horizontal_theta_mass_k_kg_m2_s,
      third.horizontal_theta_mass_k_kg_m2_s);
  budget.integrated_tracer_mass_source_kg_m2 += rk_integral(
      first.horizontal_tracer_mass_kg_m2_s, second.horizontal_tracer_mass_kg_m2_s,
      third.horizontal_tracer_mass_kg_m2_s);
}
}  // namespace

AtmosphericHybridCoordinate make_vertical_coordinate(const ExperimentConfig& config) {
  if (config.kind != ExperimentKind::kVerticalColumn) {
    throw std::invalid_argument(
        "vertical coordinate requires a vertical_column config");
  }
  return {{config.vertical.a_half_pa, config.vertical.b_half},
          config.vertical.minimum_surface_pressure_pa,
          config.vertical.maximum_surface_pressure_pa,
          config.vertical.minimum_pressure_thickness_pa};
}

VerticalColumnState make_vertical_column_initial_state(
    const ExperimentConfig& config, const AtmosphericHybridCoordinate& coordinate) {
  const auto geometry = coordinate.geometry(
      config.vertical.surface_pressure_pa, config.planet.gravity_m_s2,
      config.planet.gas_constant_j_kg_k, config.planet.heat_capacity_cp_j_kg_k,
      config.planet.reference_pressure_pa);
  VerticalColumnState state;
  state.time_s = config.run.start_time_s;
  state.surface_pressure_pa = config.vertical.surface_pressure_pa;
  state.potential_temperature_mass_k_kg_m2.resize(coordinate.levels());
  state.tracer_mass_kg_m2.resize(coordinate.levels());
  for (std::size_t k = 0; k < coordinate.levels(); ++k) {
    Real theta =
        config.vertical.test_case == VerticalTestCase::kIsothermal
            ? thermodynamics::potential_temperature_from_temperature(
                  config.vertical.initial_temperature_k, geometry.pressure_full_pa[k],
                  config.planet.gas_constant_j_kg_k,
                  config.planet.heat_capacity_cp_j_kg_k,
                  config.planet.reference_pressure_pa)
            : config.vertical.initial_potential_temperature_k;
    Real tracer = 0.1;
    if (config.vertical.test_case == VerticalTestCase::kManufacturedTransport) {
      const Real phase = 2.0 * std::numbers::pi * (static_cast<Real>(k) + 0.5) /
                         static_cast<Real>(coordinate.levels());
      theta += 5.0 * std::sin(phase);
      tracer += 0.01 * std::sin(phase);
    }
    state.potential_temperature_mass_k_kg_m2[k] = geometry.air_mass_kg_m2[k] * theta;
    state.tracer_mass_kg_m2[k] = geometry.air_mass_kg_m2[k] * tracer;
  }
  validate(config, coordinate, state);
  return state;
}

VerticalColumnResult run_vertical_column(
    const ExperimentConfig& config, std::optional<VerticalColumnState> initial_state,
    std::optional<std::uint64_t> stop_after_step) {
  config.validate();
  const auto coordinate = make_vertical_coordinate(config);
  VerticalColumnState state =
      initial_state.value_or(make_vertical_column_initial_state(config, coordinate));
  validate(config, coordinate, state);
  const auto initial_geometry = coordinate.geometry(
      state.surface_pressure_pa, config.planet.gravity_m_s2,
      config.planet.gas_constant_j_kg_k, config.planet.heat_capacity_cp_j_kg_k,
      config.planet.reference_pressure_pa);
  VerticalColumnBudget budget{
      .initial_dry_mass_kg_m2 = sum(initial_geometry.air_mass_kg_m2),
      .initial_potential_temperature_mass_k_kg_m2 =
          sum(state.potential_temperature_mass_k_kg_m2),
      .initial_tracer_mass_kg_m2 = sum(state.tracer_mass_kg_m2)};
  VerticalColumnResult result{.state = state,
                              .reached_end_time = true,
                              .maximum_cfl = 0.0,
                              .samples = {},
                              .budget = budget,
                              .mass_flux = {}};
  while (state.time_s < config.run.end_time_s) {
    if (stop_after_step.has_value() && state.step >= *stop_after_step) {
      result.reached_end_time = false;
      break;
    }
    const auto geometry = coordinate.geometry(
        state.surface_pressure_pa, config.planet.gravity_m_s2,
        config.planet.gas_constant_j_kg_k, config.planet.heat_capacity_cp_j_kg_k,
        config.planet.reference_pressure_pa);
    const auto r1 = rhs(config, coordinate, state);
    Real dt = std::min(config.run.end_time_s - state.time_s,
                       vertical_stable_time_step(
                           geometry.air_mass_kg_m2, r1.mass_flux.interface_flux_kg_m2_s,
                           r1.horizontal_air_mass_kg_m2_s, config.vertical.cfl,
                           config.run.time_step_s));
    if (config.vertical.test_case == VerticalTestCase::kIsothermal ||
        config.vertical.test_case == VerticalTestCase::kDryAdiabatic) {
      state.time_s += dt;
      ++state.step;
      if (state.step % config.diagnostics.interval_steps == 0 ||
          state.time_s == config.run.end_time_s) {
        result.samples.push_back({state.time_s, state.step, state, r1.mass_flux,
                                  result.maximum_cfl, budget});
      }
      continue;
    }
    dt = pressure_limited_time_step(config, state.surface_pressure_pa,
                                    r1.surface_pressure_pa_s, dt);
    VerticalColumnState next_state;
    VerticalRhs r2{};
    VerticalRhs r3{};
    Real step_cfl = 0.0;
    while (true) {
      auto first_stage = plus(state, r1, dt);
      first_stage.time_s = state.time_s + dt;
      if (!pressure_is_in_range(config, first_stage.surface_pressure_pa)) {
        halve_time_step(state, dt);
        continue;
      }
      validate(config, coordinate, first_stage);
      r2 = rhs(config, coordinate, first_stage);
      const auto first_stage_geometry = coordinate.geometry(
          first_stage.surface_pressure_pa, config.planet.gravity_m_s2,
          config.planet.gas_constant_j_kg_k, config.planet.heat_capacity_cp_j_kg_k,
          config.planet.reference_pressure_pa);
      const Real second_cfl = vertical_maximum_cfl(first_stage_geometry.air_mass_kg_m2,
                                                   r2.mass_flux.interface_flux_kg_m2_s,
                                                   r2.horizontal_air_mass_kg_m2_s, dt);
      if (second_cfl > config.vertical.cfl) {
        halve_time_step(state, dt);
        continue;
      }

      auto second_stage = blend(state, 0.75, plus(first_stage, r2, dt), 0.25);
      second_stage.time_s = state.time_s + 0.5 * dt;
      if (!pressure_is_in_range(config, second_stage.surface_pressure_pa)) {
        halve_time_step(state, dt);
        continue;
      }
      validate(config, coordinate, second_stage);
      r3 = rhs(config, coordinate, second_stage);
      const auto second_stage_geometry = coordinate.geometry(
          second_stage.surface_pressure_pa, config.planet.gravity_m_s2,
          config.planet.gas_constant_j_kg_k, config.planet.heat_capacity_cp_j_kg_k,
          config.planet.reference_pressure_pa);
      const Real third_cfl = vertical_maximum_cfl(second_stage_geometry.air_mass_kg_m2,
                                                  r3.mass_flux.interface_flux_kg_m2_s,
                                                  r3.horizontal_air_mass_kg_m2_s, dt);
      if (third_cfl > config.vertical.cfl) {
        halve_time_step(state, dt);
        continue;
      }

      next_state = blend(state, 1.0 / 3.0, plus(second_stage, r3, dt), 2.0 / 3.0);
      next_state.time_s = state.time_s + dt;
      next_state.step = state.step + 1;
      if (!pressure_is_in_range(config, next_state.surface_pressure_pa)) {
        halve_time_step(state, dt);
        continue;
      }
      validate(config, coordinate, next_state);
      step_cfl = std::max({vertical_maximum_cfl(geometry.air_mass_kg_m2,
                                                r1.mass_flux.interface_flux_kg_m2_s,
                                                r1.horizontal_air_mass_kg_m2_s, dt),
                           second_cfl, third_cfl});
      break;
    }
    accumulate_budget(budget, r1, r2, r3, dt);
    state = std::move(next_state);
    result.maximum_cfl = std::max(result.maximum_cfl, step_cfl);
    if (state.step % config.diagnostics.interval_steps == 0 ||
        state.time_s == config.run.end_time_s) {
      result.samples.push_back({state.time_s, state.step, state,
                                rhs(config, coordinate, state).mass_flux,
                                result.maximum_cfl, budget});
    }
  }
  result.state = state;
  result.budget = budget;
  result.mass_flux = rhs(config, coordinate, state).mass_flux;
  return result;
}

}  // namespace mps
