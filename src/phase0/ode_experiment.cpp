#include "myplanetsim/phase0/ode_experiment.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

#include "myplanetsim/core/validation.hpp"
#include "myplanetsim/numerics/explicit_steppers.hpp"

namespace mps {
namespace {

template <typename Stepper>
[[nodiscard]] OdeState integrate(const ExperimentConfig& config, OdeState state,
                                 const std::optional<std::uint64_t> stop_after) {
  Stepper stepper(1);
  std::vector<Real> values{state.value};
  const auto rhs = [&config](const Real, const std::span<const Real> current,
                             const std::span<Real> tendency) {
    tendency[0] = -config.ode.decay_rate_s_1 * current[0];
  };

  while (state.time_s < config.run.end_time_s &&
         (!stop_after.has_value() || state.step < *stop_after)) {
    const Real time_step_s =
        std::min(config.run.time_step_s, config.run.end_time_s - state.time_s);
    stepper.step(state.time_s, time_step_s, values, rhs);
    state.time_s += time_step_s;
    const Real time_tolerance =
        8.0 * std::numeric_limits<Real>::epsilon() *
        std::max({1.0, std::abs(state.time_s), std::abs(config.run.end_time_s)});
    if (config.run.end_time_s - state.time_s <= time_tolerance) {
      state.time_s = config.run.end_time_s;
    }
    ++state.step;
  }
  state.value = values[0];
  return state;
}

}  // namespace

IntegratorKind parse_integrator(const std::string_view name) {
  if (name == "euler") {
    return IntegratorKind::kForwardEuler;
  }
  if (name == "ssprk3") {
    return IntegratorKind::kSspRk3;
  }
  throw std::invalid_argument("unknown integrator: " + std::string(name));
}

std::string_view integrator_name(const IntegratorKind integrator) {
  switch (integrator) {
    case IntegratorKind::kForwardEuler:
      return "euler";
    case IntegratorKind::kSspRk3:
      return "ssprk3";
  }
  throw std::invalid_argument("invalid integrator enum value");
}

OdeResult run_ode_experiment(const ExperimentConfig& config,
                             const IntegratorKind integrator,
                             const std::optional<OdeState> initial_state,
                             const std::optional<std::uint64_t> stop_after_step) {
  config.validate();
  OdeState state = initial_state.value_or(OdeState{
      .time_s = config.run.start_time_s,
      .step = 0,
      .value = config.ode.initial_value,
  });
  require_finite(state.time_s, "ODE state time_s");
  require_finite(state.value, "ODE state value");
  if (state.time_s < config.run.start_time_s || state.time_s > config.run.end_time_s) {
    throw std::invalid_argument("ODE state time is outside the configured interval");
  }
  if (stop_after_step.has_value() && *stop_after_step < state.step) {
    throw std::invalid_argument("stop-after step precedes the initial state step");
  }

  switch (integrator) {
    case IntegratorKind::kForwardEuler:
      state = integrate<ForwardEuler>(config, state, stop_after_step);
      break;
    case IntegratorKind::kSspRk3:
      state = integrate<SspRk3>(config, state, stop_after_step);
      break;
  }

  const Real elapsed_s = state.time_s - config.run.start_time_s;
  const Real exact_value =
      config.ode.initial_value * std::exp(-config.ode.decay_rate_s_1 * elapsed_s);
  return OdeResult{
      .state = state,
      .exact_value = exact_value,
      .absolute_error = std::abs(state.value - exact_value),
      .reached_end_time = state.time_s == config.run.end_time_s,
  };
}

}  // namespace mps
