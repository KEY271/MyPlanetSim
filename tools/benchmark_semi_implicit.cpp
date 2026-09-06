#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string_view>

#if defined(__unix__) || defined(__APPLE__)
#include <sys/resource.h>
#endif

#include "myplanetsim/config/experiment_config.hpp"
#include "myplanetsim/dynamics/dry_hydrostatic_driver.hpp"

namespace {

struct CommandLine {
  std::string_view config_path;
  std::size_t requested_steps = 10;
};

[[nodiscard]] std::size_t parse_positive_count(const std::string_view text) {
  std::size_t value = 0;
  for (const char character : text) {
    if (character < '0' || character > '9')
      throw std::invalid_argument("step count must be a positive integer");
    value = value * 10 + static_cast<std::size_t>(character - '0');
  }
  if (value == 0) throw std::invalid_argument("step count must be positive");
  return value;
}

[[nodiscard]] CommandLine parse_command_line(const int argc, char** argv) {
  if (argc != 2 && argc != 4)
    throw std::invalid_argument(
        "usage: benchmark_semi_implicit CONFIG [--steps COUNT]");
  CommandLine result{.config_path = argv[1]};
  if (argc == 4) {
    if (std::string_view(argv[2]) != "--steps")
      throw std::invalid_argument(
          "usage: benchmark_semi_implicit CONFIG [--steps COUNT]");
    result.requested_steps = parse_positive_count(argv[3]);
  }
  return result;
}

[[nodiscard]] std::size_t peak_rss_bytes() noexcept {
#if defined(__unix__) || defined(__APPLE__)
  rusage usage{};
  if (getrusage(RUSAGE_SELF, &usage) != 0) return 0;
#if defined(__APPLE__)
  return static_cast<std::size_t>(usage.ru_maxrss);
#else
  return static_cast<std::size_t>(usage.ru_maxrss) * 1024;
#endif
#else
  return 0;
#endif
}

}  // namespace

int main(const int argc, char** argv) {
  try {
    const auto command_line = parse_command_line(argc, argv);
    const auto config = mps::load_experiment_config(command_line.config_path);
    if (config.kind != mps::ExperimentKind::kDryHydrostatic ||
        config.dry_hydrostatic.time_integrator !=
            mps::DryHydrostaticTimeIntegrator::kSemiImplicit)
      throw std::invalid_argument(
          "benchmark requires a semi-implicit dry-hydrostatic configuration");
    const mps::DryHydrostaticDriver driver(config);
    auto state = driver.initial_state();
    const auto initial_step = state.step;
    const auto initial_time_s = state.time_s;
    const auto end_time_s =
        state.time_s +
        static_cast<mps::Real>(command_line.requested_steps) * config.run.time_step_s;

    std::size_t diagnostic_steps = 0;
    std::size_t linear_iterations_total = 0;
    std::size_t linear_iterations_maximum = 0;
    std::size_t retries = 0;
    mps::Real accepted_dt_sum = 0.0;
    mps::Real accepted_dt_minimum = std::numeric_limits<mps::Real>::infinity();
    mps::Real accepted_dt_maximum = 0.0;
    mps::Real rhs_wall_seconds = 0.0;
    mps::Real linear_wall_seconds = 0.0;
    const auto start = std::chrono::steady_clock::now();
    driver.advance(
        state, end_time_s,
        [&](const mps::DryHydrostaticState& sampled, const mps::DryHydrostaticDerived*,
            const mps::DryHydrostaticStepDiagnostics& diagnostics) {
          if (sampled.step == initial_step) return;
          ++diagnostic_steps;
          linear_iterations_total += diagnostics.linear_iterations_total;
          linear_iterations_maximum = std::max(linear_iterations_maximum,
                                               diagnostics.linear_iterations_maximum);
          retries += diagnostics.retry_count;
          accepted_dt_sum += diagnostics.accepted_time_step_s;
          accepted_dt_minimum =
              std::min(accepted_dt_minimum, diagnostics.accepted_time_step_s);
          accepted_dt_maximum =
              std::max(accepted_dt_maximum, diagnostics.accepted_time_step_s);
          rhs_wall_seconds += diagnostics.wall_seconds_rhs;
          linear_wall_seconds += diagnostics.wall_seconds_linear_solve;
        });
    const mps::Real elapsed_s =
        std::chrono::duration<mps::Real>(std::chrono::steady_clock::now() - start)
            .count();
    const auto completed_steps = state.step - initial_step;
    const mps::Real advanced_model_s = state.time_s - initial_time_s;
    if (completed_steps == 0 || diagnostic_steps != completed_steps ||
        !(advanced_model_s > 0.0))
      throw std::runtime_error("benchmark did not complete an accepted step");

    std::cout << std::setprecision(10) << "benchmark = semi_implicit_dry_hydrostatic\n"
              << "requested_steps = " << command_line.requested_steps << '\n'
              << "accepted_steps = " << completed_steps << '\n'
              << "advanced_model_s = " << advanced_model_s << '\n'
              << "elapsed_s = " << elapsed_s << '\n'
              << "seconds_per_model_day = " << elapsed_s * 86400.0 / advanced_model_s
              << '\n'
              << "accepted_dt_mean_s = "
              << accepted_dt_sum / static_cast<mps::Real>(completed_steps) << '\n'
              << "accepted_dt_minimum_s = " << accepted_dt_minimum << '\n'
              << "accepted_dt_maximum_s = " << accepted_dt_maximum << '\n'
              << "linear_iterations_total = " << linear_iterations_total << '\n'
              << "linear_iterations_maximum = " << linear_iterations_maximum << '\n'
              << "retry_count = " << retries << '\n'
              << "wall_seconds_rhs = " << rhs_wall_seconds << '\n'
              << "wall_seconds_linear_solve = " << linear_wall_seconds << '\n'
              << "wall_seconds_unattributed = "
              << std::max(0.0, elapsed_s - rhs_wall_seconds - linear_wall_seconds)
              << '\n'
              << "peak_rss_bytes = " << peak_rss_bytes() << '\n';
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "semi-implicit benchmark error: " << error.what() << '\n';
    return 1;
  }
}
