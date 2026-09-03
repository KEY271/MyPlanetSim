#include <charconv>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>

#include "myplanetsim/config/experiment_config.hpp"
#include "myplanetsim/dynamics/shallow_water_driver.hpp"
#include "myplanetsim/grid/cubed_sphere_grid.hpp"
#include "myplanetsim/io/checkpoint.hpp"
#include "myplanetsim/io/run_metadata.hpp"
#include "myplanetsim/phase0/ode_experiment.hpp"
#include "myplanetsim/transport/spherical_transport.hpp"

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

void write_transport_result(std::ostream& output, const mps::TransportResult& result) {
  const auto& diagnostics = result.diagnostics;
  output << std::setprecision(std::numeric_limits<mps::Real>::max_digits10)
         << "result.integrator = ssprk3\n"
         << "result.status = " << (result.reached_end_time ? "complete" : "stopped")
         << '\n'
         << "result.time_s = " << result.state.time_s << '\n'
         << "result.step = " << result.state.step << '\n'
         << "diagnostics.initial_mass = " << diagnostics.initial_mass << '\n'
         << "diagnostics.final_mass = " << diagnostics.final_mass << '\n'
         << "diagnostics.relative_mass_drift = " << diagnostics.relative_mass_drift
         << '\n'
         << "diagnostics.l1_error = " << diagnostics.l1_error << '\n'
         << "diagnostics.l2_error = " << diagnostics.l2_error << '\n'
         << "diagnostics.linf_error = " << diagnostics.linf_error << '\n'
         << "diagnostics.minimum = " << diagnostics.minimum << '\n'
         << "diagnostics.maximum = " << diagnostics.maximum << '\n'
         << "diagnostics.seam_rms_error = " << diagnostics.seam_rms_error << '\n'
         << "diagnostics.corner_rms_error = " << diagnostics.corner_rms_error << '\n'
         << "diagnostics.filament_preservation = " << diagnostics.filament_preservation
         << '\n'
         << "diagnostics.unmixing = " << diagnostics.unmixing << '\n'
         << "diagnostics.overshooting = " << diagnostics.overshooting << '\n'
         << "diagnostics.real_mixing = " << diagnostics.real_mixing << '\n'
         << "diagnostics.limiter_activations = " << diagnostics.limiter_activations
         << '\n';
}

void write_transport_snapshot(const mps::ExperimentConfig& config,
                              const mps::TransportResult& result) {
  const std::filesystem::path directory(config.output_directory);
  std::filesystem::create_directories(directory);
  std::ofstream output(directory / "tracer.csv", std::ios::trunc);
  if (!output) {
    throw std::runtime_error("unable to open transport CSV snapshot");
  }
  output << "panel,i,j,tracer\n";
  const mps::CubedSphereGrid grid(config.grid.cells_per_panel, config.planet.radius_m);
  output << std::setprecision(std::numeric_limits<mps::Real>::max_digits10);
  for (std::size_t index = 0; index < grid.cell_count(); ++index) {
    const auto cell = grid.cell_id(index);
    output << mps::panel_name(cell.panel) << ',' << cell.i << ',' << cell.j << ','
           << result.state.tracer[index] << '\n';
  }
}

void write_shallow_water_result(std::ostream& output,
                                const mps::ShallowWaterResult& result) {
  output << std::setprecision(std::numeric_limits<mps::Real>::max_digits10)
         << "result.integrator = ssprk3\n"
         << "result.status = " << (result.reached_end_time ? "complete" : "stopped")
         << '\n'
         << "result.time_s = " << result.state.time_s << '\n'
         << "result.step = " << result.state.step << '\n'
         << "diagnostics.mass = " << result.final_diagnostics.mass << '\n'
         << "diagnostics.energy = " << result.final_diagnostics.energy << '\n'
         << "diagnostics.potential_enstrophy = "
         << result.final_diagnostics.potential_enstrophy << '\n'
         << "diagnostics.axial_angular_momentum = "
         << result.final_diagnostics.axial_angular_momentum << '\n'
         << "diagnostics.minimum_depth = " << result.final_diagnostics.minimum_depth
         << '\n'
         << "diagnostics.maximum_depth = " << result.final_diagnostics.maximum_depth
         << '\n'
         << "diagnostics.maximum_cfl = " << result.maximum_cfl << '\n';
}

void write_shallow_water_snapshot(const mps::ExperimentConfig& config,
                                  const mps::ShallowWaterResult& result) {
  const std::filesystem::path directory(config.output_directory);
  std::filesystem::create_directories(directory);
  std::ofstream output(directory / "shallow_water.csv", std::ios::trunc);
  if (!output) {
    throw std::runtime_error("unable to open shallow-water CSV snapshot");
  }
  output << "panel,i,j,depth_m,momentum_x,momentum_y,momentum_z\n";
  const mps::CubedSphereGrid grid(config.grid.cells_per_panel, config.planet.radius_m);
  output << std::setprecision(std::numeric_limits<mps::Real>::max_digits10);
  for (std::size_t index = 0; index < grid.cell_count(); ++index) {
    const auto cell = grid.cell_id(index);
    const auto momentum = result.state.momentum[index];
    output << mps::panel_name(cell.panel) << ',' << cell.i << ',' << cell.j << ','
           << result.state.depth[index] << ',' << momentum.x << ',' << momentum.y << ','
           << momentum.z << '\n';
  }
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

    if (config.kind == mps::ExperimentKind::kOde) {
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
    } else if (config.kind == mps::ExperimentKind::kSphereTransport) {
      if (command_line.integrator != mps::IntegratorKind::kSspRk3) {
        throw std::invalid_argument("sphere transport supports only ssprk3");
      }
      std::optional<mps::TransportState> initial_state;
      if (command_line.restart_path.has_value()) {
        auto checkpoint =
            mps::read_checkpoint_file(*command_line.restart_path, fingerprint);
        initial_state = mps::TransportState{checkpoint.time_s, checkpoint.step,
                                            std::move(checkpoint.state)};
      }
      const auto result = mps::run_spherical_transport(config, std::move(initial_state),
                                                       command_line.stop_after_step);
      if (command_line.checkpoint_path.has_value()) {
        mps::write_checkpoint_file(
            *command_line.checkpoint_path,
            mps::Checkpoint{result.state.time_s, result.state.step, result.state.tracer,
                            fingerprint});
      }
      write_transport_snapshot(config, result);
      mps::write_run_metadata(std::cout, mps::make_run_metadata(config), config);
      write_transport_result(std::cout, result);
    } else {
      if (command_line.integrator != mps::IntegratorKind::kSspRk3) {
        throw std::invalid_argument("shallow water supports only ssprk3");
      }
      const mps::CubedSphereGrid grid(config.grid.cells_per_panel,
                                      config.planet.radius_m);
      std::optional<mps::ShallowWaterState> initial_state;
      if (command_line.restart_path.has_value()) {
        auto checkpoint = mps::read_checkpoint_file(
            *command_line.restart_path, fingerprint, mps::kShallowWaterCheckpointLayout,
            4 * grid.cell_count());
        initial_state = mps::unflatten_shallow_water_state(
            checkpoint.time_s, checkpoint.step, checkpoint.state, grid.cell_count());
      }
      const auto result = mps::run_shallow_water(config, std::move(initial_state),
                                                 command_line.stop_after_step);
      if (command_line.checkpoint_path.has_value()) {
        mps::write_checkpoint_file(
            *command_line.checkpoint_path,
            mps::Checkpoint{
                .time_s = result.state.time_s,
                .step = result.state.step,
                .state = mps::flatten_shallow_water_state(result.state),
                .config_fingerprint = fingerprint,
                .layout_id = std::string(mps::kShallowWaterCheckpointLayout),
            });
      }
      write_shallow_water_snapshot(config, result);
      mps::write_run_metadata(std::cout, mps::make_run_metadata(config), config);
      write_shallow_water_result(std::cout, result);
    }
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << '\n';
    return 1;
  }
}
