#pragma once

#include <cstdint>
#include <optional>
#include <string_view>

#include "myplanetsim/config/experiment_config.hpp"
#include "myplanetsim/core/types.hpp"

namespace mps {

enum class IntegratorKind { kForwardEuler, kSspRk3 };

struct OdeState {
  Real time_s;
  std::uint64_t step;
  Real value;
};

struct OdeResult {
  OdeState state;
  Real exact_value;
  Real absolute_error;
  bool reached_end_time;
};

[[nodiscard]] IntegratorKind parse_integrator(std::string_view name);
[[nodiscard]] std::string_view integrator_name(IntegratorKind integrator);
[[nodiscard]] OdeResult run_ode_experiment(
    const ExperimentConfig& config, IntegratorKind integrator,
    std::optional<OdeState> initial_state = std::nullopt,
    std::optional<std::uint64_t> stop_after_step = std::nullopt);

}  // namespace mps
