#include <charconv>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>

#include "myplanetsim/config/experiment_config.hpp"
#include "myplanetsim/io/checkpoint.hpp"
#include "myplanetsim/io/run_metadata.hpp"
#include "myplanetsim/phase0/ode_experiment.hpp"

namespace {

struct CommandLine {
  std::filesystem::path config_path;
  mps::IntegratorKind integrator = mps::IntegratorKind::kSspRk3;
  std::optional<std::filesystem::path> checkpoint_path;
  std::optional<std::filesystem::path> restart_path;
  std::optional<std::uint64_t> stop_after_step;
};

void print_usage(std::ostream& output) {
  output << "Usage: my_planet_sim --config PATH [options]\n"
         << "Options:\n"
         << "  --integrator euler|ssprk3  Time integrator (default: ssprk3)\n"
         << "  --checkpoint PATH          Write the final or stopped state\n"
         << "  --restart PATH             Continue from a checkpoint\n"
         << "  --stop-after-step N        Stop after absolute step N\n"
         << "  --help                     Show this help\n";
}

[[nodiscard]] std::uint64_t parse_step(const std::string_view text) {
  std::uint64_t value = 0;
  const auto result = std::from_chars(text.data(), text.data() + text.size(), value);
  if (result.ec != std::errc{} || result.ptr != text.data() + text.size()) {
    throw std::invalid_argument("invalid --stop-after-step value");
  }
  return value;
}

[[nodiscard]] CommandLine parse_command_line(const int argc, const char* const argv[]) {
  CommandLine command_line;
  for (int index = 1; index < argc; ++index) {
    const std::string_view argument(argv[index]);
    if (argument == "--help") {
      print_usage(std::cout);
      throw std::runtime_error("help requested");
    }
    if (index + 1 >= argc) {
      throw std::invalid_argument("missing value after " + std::string(argument));
    }
    const std::string_view value(argv[++index]);
    if (argument == "--config") {
      command_line.config_path = value;
    } else if (argument == "--integrator") {
      command_line.integrator = mps::parse_integrator(value);
    } else if (argument == "--checkpoint") {
      command_line.checkpoint_path = value;
    } else if (argument == "--restart") {
      command_line.restart_path = value;
    } else if (argument == "--stop-after-step") {
      command_line.stop_after_step = parse_step(value);
    } else {
      throw std::invalid_argument("unknown option: " + std::string(argument));
    }
  }
  if (command_line.config_path.empty()) {
    throw std::invalid_argument("--config is required");
  }
  return command_line;
}

void write_result(std::ostream& output, const mps::OdeResult& result,
                  const mps::IntegratorKind integrator) {
  output << std::setprecision(std::numeric_limits<mps::Real>::max_digits10)
         << "result.integrator = " << mps::integrator_name(integrator) << '\n'
         << "result.status = " << (result.reached_end_time ? "complete" : "stopped")
         << '\n'
         << "result.time_s = " << result.state.time_s << '\n'
         << "result.step = " << result.state.step << '\n'
         << "result.value = " << result.state.value << '\n'
         << "result.exact_value = " << result.exact_value << '\n'
         << "result.absolute_error = " << result.absolute_error << '\n';
}

}  // namespace

int main(const int argc, const char* const argv[]) {
  if (argc == 2 && std::string_view(argv[1]) == "--help") {
    print_usage(std::cout);
    return 0;
  }

  try {
    const auto command_line = parse_command_line(argc, argv);
    const auto config = mps::load_experiment_config(command_line.config_path);
    const auto fingerprint = mps::config_fingerprint(config);

    std::optional<mps::OdeState> initial_state;
    if (command_line.restart_path.has_value()) {
      const auto checkpoint =
          mps::read_checkpoint_file(*command_line.restart_path, fingerprint);
      if (checkpoint.state.size() != 1) {
        throw std::runtime_error("ODE checkpoint must contain exactly one value");
      }
      initial_state = mps::OdeState{.time_s = checkpoint.time_s,
                                    .step = checkpoint.step,
                                    .value = checkpoint.state[0]};
    }

    const auto result = mps::run_ode_experiment(
        config, command_line.integrator, initial_state, command_line.stop_after_step);
    if (command_line.checkpoint_path.has_value()) {
      mps::write_checkpoint_file(*command_line.checkpoint_path,
                                 mps::Checkpoint{.time_s = result.state.time_s,
                                                 .step = result.state.step,
                                                 .state = {result.state.value},
                                                 .config_fingerprint = fingerprint});
    }

    mps::write_run_metadata(std::cout, mps::make_run_metadata(config), config);
    write_result(std::cout, result, command_line.integrator);
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << '\n';
    return 1;
  }
}
