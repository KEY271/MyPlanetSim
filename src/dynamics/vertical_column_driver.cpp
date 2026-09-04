#include "myplanetsim/dynamics/vertical_column_driver.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <stdexcept>

#include "myplanetsim/thermodynamics/dry_thermodynamics.hpp"

namespace mps {
namespace {
struct VerticalRhs {
  Real surface_pressure_pa_s;
  std::vector<Real> horizontal_air_mass_kg_m2_s;
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
    const Real phase = 2.0 * std::numbers::pi * state.time_s /
                       std::max(config.run.end_time_s - config.run.start_time_s, 1.0);
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
  const auto geometry = coordinate.geometry(
      state.surface_pressure_pa, config.planet.gravity_m_s2,
      config.planet.gas_constant_j_kg_k, config.planet.heat_capacity_cp_j_kg_k,
      config.planet.reference_pressure_pa);
  validate_vertical_column_state(state, geometry.air_mass_kg_m2,
                                 config.vertical.temperature_floor_k,
                                 geometry.exner_full);
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
  VerticalColumnResult result{
      .state = state, .reached_end_time = true, .maximum_cfl = 0.0};
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
    const Real dt =
        std::min(config.run.end_time_s - state.time_s,
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
        result.samples.push_back(
            {state.time_s, state.step, state, r1.mass_flux, result.maximum_cfl});
      }
      continue;
    }
    auto stage = plus(state, r1, dt);
    stage.time_s = state.time_s + dt;
    validate(config, coordinate, stage);
    const auto r2 = rhs(config, coordinate, stage);
    stage = blend(state, 0.75, plus(stage, r2, dt), 0.25);
    stage.time_s = state.time_s + 0.5 * dt;
    validate(config, coordinate, stage);
    const auto r3 = rhs(config, coordinate, stage);
    state = blend(state, 1.0 / 3.0, plus(stage, r3, dt), 2.0 / 3.0);
    state.time_s += dt;
    ++state.step;
    validate(config, coordinate, state);
    result.maximum_cfl =
        std::max(result.maximum_cfl, dt / config.run.time_step_s * config.vertical.cfl);
    if (state.step % config.diagnostics.interval_steps == 0 ||
        state.time_s == config.run.end_time_s) {
      result.samples.push_back({state.time_s, state.step, state,
                                rhs(config, coordinate, state).mass_flux,
                                result.maximum_cfl});
    }
  }
  result.state = state;
  return result;
}

}  // namespace mps
