#include <atomic>
#include <charconv>
#include <chrono>
#include <csignal>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "myplanetsim/config/experiment_config.hpp"
#include "myplanetsim/control/control_request.hpp"
#include "myplanetsim/diagnostics/reductions.hpp"
#include "myplanetsim/diagnostics/vertical_column_diagnostics.hpp"
#include "myplanetsim/dynamics/shallow_water_benchmarks.hpp"
#include "myplanetsim/dynamics/shallow_water_driver.hpp"
#include "myplanetsim/dynamics/vertical_column_driver.hpp"
#include "myplanetsim/grid/cubed_sphere_grid.hpp"
#include "myplanetsim/io/checkpoint.hpp"
#include "myplanetsim/io/frame.hpp"
#include "myplanetsim/io/run_metadata.hpp"
#include "myplanetsim/phase0/ode_experiment.hpp"
#include "myplanetsim/transport/spherical_transport.hpp"
#include "myplanetsim/vertical/hydrostatic_column.hpp"

namespace {

struct CommandLine {
  std::filesystem::path config_path;
  mps::IntegratorKind integrator = mps::IntegratorKind::kSspRk3;
  std::optional<std::filesystem::path> checkpoint_path;
  std::optional<std::filesystem::path> restart_path;
  std::optional<std::uint64_t> stop_after_step;
  std::optional<std::filesystem::path> control_request_path;
  bool event_stream_ndjson = false;
  bool describe_control = false;
};

void print_usage(std::ostream& output) {
  output << "Usage: my_planet_sim --config PATH [options]\n"
         << "Options:\n"
         << "  --integrator euler|ssprk3  Time integrator (default: ssprk3)\n"
         << "  --checkpoint PATH          Write the final or stopped state\n"
         << "  --restart PATH             Continue from a checkpoint\n"
         << "  --stop-after-step N        Stop after absolute step N\n"
         << "  --control-request PATH     Run a validated machine control request\n"
         << "  --event-stream ndjson      Emit machine-readable lifecycle events\n"
         << "  --describe-control         Describe the machine control protocol\n"
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
    if (argument == "--describe-control") {
      command_line.describe_control = true;
      continue;
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
    } else if (argument == "--control-request") {
      command_line.control_request_path = value;
    } else if (argument == "--event-stream") {
      if (value != "ndjson") {
        throw std::invalid_argument("unsupported --event-stream format");
      }
      command_line.event_stream_ndjson = true;
    } else {
      throw std::invalid_argument("unknown option: " + std::string(argument));
    }
  }
  if (command_line.describe_control) {
    return command_line;
  }
  if (command_line.config_path.empty()) {
    throw std::invalid_argument("--config is required");
  }
  return command_line;
}

std::atomic_bool g_cancel_requested = false;

void request_cancellation(const int) noexcept { g_cancel_requested.store(true); }

bool cancellation_requested(void*) noexcept { return g_cancel_requested.load(); }

[[nodiscard]] std::string json_escape(const std::string_view value) {
  std::string escaped;
  escaped.reserve(value.size());
  for (const char character : value) {
    if (character == '\\' || character == '"') {
      escaped += '\\';
      escaped += character;
    } else if (character == '\n') {
      escaped += "\\n";
    } else if (character == '\r') {
      escaped += "\\r";
    } else {
      escaped += character;
    }
  }
  return escaped;
}

struct MachineEventWriter {
  std::string run_id;
  std::uint64_t sequence = 0;

  void emit(const std::string_view type, const std::string_view fields = {}) {
    std::cout << "{\"protocolVersion\":1,\"runId\":\"" << json_escape(run_id)
              << "\",\"sequence\":" << sequence++ << ",\"type\":\"" << type << '"';
    if (!fields.empty()) {
      std::cout << ',' << fields;
    }
    std::cout << "}\n" << std::flush;
    if (!std::cout) {
      throw std::runtime_error("failed while writing machine event stream");
    }
  }
};

struct MachineFrameContext {
  std::filesystem::path directory;
  std::string config_fingerprint;
  mps::Index cells_per_panel = 0;
  MachineEventWriter* events = nullptr;
  std::uint64_t frame_sequence = 0;
};

void write_machine_frame(const mps::ShallowWaterState& state, void* context) {
  auto& frame = *static_cast<MachineFrameContext*>(context);
  const auto sequence = frame.frame_sequence++;
  const std::string filename = "frame_" + std::to_string(sequence) + ".bin";
  const auto path = frame.directory / filename;
  mps::write_frame_file(path,
                        {.cells_per_panel = frame.cells_per_panel,
                         .time_s = state.time_s,
                         .step = state.step,
                         .config_fingerprint = frame.config_fingerprint},
                        state);
  frame.events->emit(
      "frame.ready",
      "\"frameSequence\":" + std::to_string(sequence) + ",\"timeSeconds\":" +
          std::to_string(state.time_s) + ",\"step\":" + std::to_string(state.step) +
          ",\"relativePath\":\"" + json_escape(filename) +
          "\",\"byteLength\":" + std::to_string(std::filesystem::file_size(path)));
}

void print_control_description() {
  std::cout << "{\"protocolVersion\":1,\"supportedEdits\":[\"gaussian_depth\"],"
               "\"maxEditCount\":64,\"maxEndTimeSeconds\":31536000,"
               "\"maxTimeStepSeconds\":86400,\"maxFrameIntervalSteps\":1000000}\n";
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
                                const mps::ShallowWaterResult& result,
                                const mps::Real wall_time_s,
                                const std::size_t cell_count) {
  const auto& initial = result.initial_diagnostics;
  const auto& final_state = result.final_diagnostics;
  const auto drift = [](const mps::Real final_value, const mps::Real initial_value) {
    return mps::diagnostics::relative_drift(final_value, initial_value, initial_value);
  };
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
         << "diagnostics.maximum_cfl = " << result.maximum_cfl << '\n'
         << "diagnostics.mass_drift = " << drift(final_state.mass, initial.mass) << '\n'
         << "diagnostics.energy_drift = " << drift(final_state.energy, initial.energy)
         << '\n'
         << "diagnostics.potential_enstrophy_drift = "
         << drift(final_state.potential_enstrophy, initial.potential_enstrophy) << '\n'
         << "diagnostics.axial_angular_momentum_drift = "
         << drift(final_state.axial_angular_momentum, initial.axial_angular_momentum)
         << '\n'
         << "cost.wall_time_s = " << wall_time_s << '\n'
         << "cost.cell_steps_per_s = "
         << static_cast<mps::Real>(cell_count) *
                static_cast<mps::Real>(result.state.step) / wall_time_s
         << '\n';
}

// Error norms against the analytic solution, where the benchmark has one.
void write_shallow_water_errors(std::ostream& output,
                                const mps::ExperimentConfig& config,
                                const mps::CubedSphereGrid& grid,
                                const mps::ShallowWaterResult& result) {
  if (config.shallow_water.test_case != mps::ShallowWaterTestCase::kWilliamson2) {
    return;
  }
  // Williamson 2 is steady, so its initial state is also the exact final state.
  const auto exact = mps::make_shallow_water_initial_state(grid, config);
  std::vector<mps::Real> weights(grid.cell_count());
  std::vector<mps::Real> velocity_error(grid.cell_count());
  const std::vector<mps::Real> zero(grid.cell_count());
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
    weights[cell] = grid.cells()[cell].area_m2;
    velocity_error[cell] =
        mps::norm(result.state.velocity(cell) - exact.velocity(cell));
  }
  const auto depth =
      mps::diagnostics::weighted_error_norms(result.state.depth, exact.depth, weights);
  const auto velocity =
      mps::diagnostics::weighted_error_norms(velocity_error, zero, weights);
  output << "error.depth_l1 = " << depth.l1 << '\n'
         << "error.depth_l2 = " << depth.l2 << '\n'
         << "error.depth_linf = " << depth.linf << '\n'
         << "error.velocity_l1 = " << velocity.l1 << '\n'
         << "error.velocity_l2 = " << velocity.l2 << '\n'
         << "error.velocity_linf = " << velocity.linf << '\n';
}

void write_shallow_water_diagnostics(const mps::ExperimentConfig& config,
                                     const mps::ShallowWaterResult& result) {
  const std::filesystem::path directory(config.output_directory);
  std::filesystem::create_directories(directory);
  std::ofstream output(directory / "diagnostics.csv", std::ios::trunc);
  if (!output) {
    throw std::runtime_error("unable to open shallow-water diagnostics CSV");
  }
  output << "time_s,step,mass,energy,potential_enstrophy,axial_angular_momentum,"
            "minimum_depth,maximum_depth,minimum_pv,maximum_pv\n";
  output << std::setprecision(std::numeric_limits<mps::Real>::max_digits10);
  for (const auto& sample : result.samples) {
    const auto& invariants = sample.invariants;
    output << sample.time_s << ',' << sample.step << ',' << invariants.mass << ','
           << invariants.energy << ',' << invariants.potential_enstrophy << ','
           << invariants.axial_angular_momentum << ',' << invariants.minimum_depth
           << ',' << invariants.maximum_depth << ',' << invariants.minimum_pv << ','
           << invariants.maximum_pv << '\n';
  }
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

void write_vertical_column_result(
    std::ostream& output, const mps::VerticalColumnResult& result,
    const mps::diagnostics::VerticalColumnDiagnostics& diagnostics) {
  output << std::setprecision(std::numeric_limits<mps::Real>::max_digits10)
         << "result.integrator = ssprk3\n"
         << "result.status = " << (result.reached_end_time ? "complete" : "stopped")
         << '\n'
         << "result.time_s = " << result.state.time_s << '\n'
         << "result.step = " << result.state.step << '\n'
         << "diagnostics.dry_mass_kg_m2 = " << diagnostics.dry_mass_kg_m2 << '\n'
         << "diagnostics.mass_residual_kg_m2 = "
         << diagnostics.dry_mass_structural_residual_kg_m2 << '\n'
         << "diagnostics.dry_mass_budget_residual_kg_m2 = "
         << diagnostics.dry_mass_budget_residual_kg_m2 << '\n'
         << "diagnostics.theta_mass_budget_residual_k_kg_m2 = "
         << diagnostics.potential_temperature_mass_budget_residual_k_kg_m2 << '\n'
         << "diagnostics.tracer_mass_budget_residual_kg_m2 = "
         << diagnostics.tracer_mass_budget_residual_kg_m2 << '\n'
         << "diagnostics.top_mass_flux_kg_m2_s = " << diagnostics.top_mass_flux_kg_m2_s
         << '\n'
         << "diagnostics.surface_mass_flux_kg_m2_s = "
         << diagnostics.surface_mass_flux_kg_m2_s << '\n'
         << "diagnostics.continuity_residual_pa_s = "
         << diagnostics.continuity_residual_pa_s << '\n'
         << "diagnostics.maximum_cfl = " << diagnostics.maximum_cfl << '\n'
         << "diagnostics.hydrostatic_l1_residual_m2_s2 = "
         << diagnostics.hydrostatic_l1_residual_m2_s2 << '\n'
         << "diagnostics.hydrostatic_l2_residual_m2_s2 = "
         << diagnostics.hydrostatic_l2_residual_m2_s2 << '\n'
         << "diagnostics.hydrostatic_linf_residual_m2_s2 = "
         << diagnostics.hydrostatic_linf_residual_m2_s2 << '\n'
         << "diagnostics.non_finite_count = " << diagnostics.non_finite_count << '\n';
}

void write_vertical_column_files(const mps::ExperimentConfig& config,
                                 const mps::AtmosphericHybridCoordinate& coordinate,
                                 const mps::VerticalColumnResult& result) {
  const std::filesystem::path directory(config.output_directory);
  std::filesystem::create_directories(directory);
  const auto geometry = coordinate.geometry(
      result.state.surface_pressure_pa, config.planet.gravity_m_s2,
      config.planet.gas_constant_j_kg_k, config.planet.heat_capacity_cp_j_kg_k,
      config.planet.reference_pressure_pa);
  std::vector<mps::Real> theta(coordinate.levels());
  for (std::size_t k = 0; k < theta.size(); ++k) {
    theta[k] =
        result.state.potential_temperature_mass_k_kg_m2[k] / geometry.air_mass_kg_m2[k];
  }
  const auto hydrostatic = mps::integrate_hydrostatic_column(
      geometry, theta, config.planet.heat_capacity_cp_j_kg_k,
      config.planet.gravity_m_s2, config.vertical.surface_geopotential_m2_s2);
  std::ofstream profile(directory / "column_profile.csv", std::ios::trunc);
  if (!profile) throw std::runtime_error("unable to open column profile CSV");
  profile << "location,index,pressure_pa,geopotential_m2_s2,height_m,air_mass_kg_m2,"
             "theta_k,temperature_k,tracer\n";
  profile << std::setprecision(std::numeric_limits<mps::Real>::max_digits10);
  for (std::size_t k = 0; k <= coordinate.levels(); ++k) {
    profile << "interface," << k << ',' << geometry.pressure_half_pa[k] << ','
            << hydrostatic.geopotential_half_m2_s2[k] << ','
            << hydrostatic.height_half_m[k] << ",,,,\n";
  }
  for (std::size_t k = 0; k < coordinate.levels(); ++k) {
    const auto tracer = result.state.tracer_mass_kg_m2[k] / geometry.air_mass_kg_m2[k];
    profile << "full," << k << ',' << geometry.pressure_full_pa[k] << ','
            << hydrostatic.geopotential_full_m2_s2[k] << ','
            << hydrostatic.height_full_m[k] << ',' << geometry.air_mass_kg_m2[k] << ','
            << theta[k] << ',' << theta[k] * geometry.exner_full[k] << ',' << tracer
            << '\n';
  }
  std::ofstream diagnostics(directory / "column_diagnostics.csv", std::ios::trunc);
  if (!diagnostics) throw std::runtime_error("unable to open column diagnostics CSV");
  diagnostics << "time_s,step,dry_mass_kg_m2,expected_dry_mass_kg_m2,"
                 "dry_mass_structural_residual_kg_m2,dry_mass_budget_residual_kg_m2,"
                 "theta_mass_k_kg_m2,theta_mass_budget_residual_k_kg_m2,"
                 "tracer_mass_kg_m2,tracer_mass_budget_residual_kg_m2,"
                 "minimum_delta_pressure_pa,maximum_delta_pressure_pa,"
                 "minimum_air_mass_kg_m2,maximum_air_mass_kg_m2,minimum_theta_k,"
                 "maximum_theta_k,minimum_temperature_k,maximum_temperature_k,"
                 "minimum_tracer,maximum_tracer,top_mass_flux_kg_m2_s,"
                 "surface_mass_flux_kg_m2_s,continuity_residual_pa_s,maximum_cfl,"
                 "hydrostatic_l1_residual_m2_s2,hydrostatic_l2_residual_m2_s2,"
                 "hydrostatic_linf_residual_m2_s2,non_finite_count\n";
  for (const auto& sample : result.samples) {
    const auto values = mps::diagnostics::diagnose_vertical_column(
        config, coordinate, sample.state, sample.mass_flux, sample.maximum_cfl,
        sample.budget);
    diagnostics << std::setprecision(std::numeric_limits<mps::Real>::max_digits10)
                << sample.time_s << ',' << sample.step << ',' << values.dry_mass_kg_m2
                << ',' << values.expected_dry_mass_kg_m2 << ','
                << values.dry_mass_structural_residual_kg_m2 << ','
                << values.dry_mass_budget_residual_kg_m2 << ','
                << values.potential_temperature_mass_k_kg_m2 << ','
                << values.potential_temperature_mass_budget_residual_k_kg_m2 << ','
                << values.tracer_mass_kg_m2 << ','
                << values.tracer_mass_budget_residual_kg_m2 << ','
                << values.minimum_delta_pressure_pa << ','
                << values.maximum_delta_pressure_pa << ','
                << values.minimum_air_mass_kg_m2 << ',' << values.maximum_air_mass_kg_m2
                << ',' << values.minimum_theta_k << ',' << values.maximum_theta_k << ','
                << values.minimum_temperature_k << ',' << values.maximum_temperature_k
                << ',' << values.minimum_tracer << ',' << values.maximum_tracer << ','
                << values.top_mass_flux_kg_m2_s << ','
                << values.surface_mass_flux_kg_m2_s << ','
                << values.continuity_residual_pa_s << ',' << values.maximum_cfl << ','
                << values.hydrostatic_l1_residual_m2_s2 << ','
                << values.hydrostatic_l2_residual_m2_s2 << ','
                << values.hydrostatic_linf_residual_m2_s2 << ','
                << values.non_finite_count << '\n';
  }
}

}  // namespace

int main(const int argc, const char* const argv[]) {
  if (argc == 2 && std::string_view(argv[1]) == "--help") {
    print_usage(std::cout);
    return 0;
  }

  std::unique_ptr<MachineEventWriter> machine_events;
  try {
    const auto command_line = parse_command_line(argc, argv);
    if (command_line.describe_control) {
      print_control_description();
      return 0;
    }
    if (command_line.event_stream_ndjson &&
        !command_line.control_request_path.has_value()) {
      throw std::invalid_argument("--event-stream requires --control-request");
    }
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
    } else if (config.kind == mps::ExperimentKind::kVerticalColumn) {
      if (command_line.integrator != mps::IntegratorKind::kSspRk3 ||
          command_line.control_request_path.has_value() ||
          command_line.event_stream_ndjson) {
        throw std::invalid_argument(
            "vertical column supports only standalone ssprk3 runs");
      }
      const auto coordinate = mps::make_vertical_coordinate(config);
      std::optional<mps::VerticalColumnState> initial_state;
      if (command_line.restart_path.has_value()) {
        const auto checkpoint = mps::read_checkpoint_file(
            *command_line.restart_path, fingerprint,
            mps::kVerticalColumnCheckpointLayout, 1 + 2 * coordinate.levels());
        initial_state = mps::unflatten_vertical_column_state(
            checkpoint.time_s, checkpoint.step, checkpoint.state, coordinate.levels());
      }
      const auto result = mps::run_vertical_column(config, std::move(initial_state),
                                                   command_line.stop_after_step);
      if (command_line.checkpoint_path.has_value()) {
        mps::write_checkpoint_file(
            *command_line.checkpoint_path,
            {.time_s = result.state.time_s,
             .step = result.state.step,
             .state = mps::flatten_vertical_column_state(result.state),
             .config_fingerprint = fingerprint,
             .layout_id = std::string(mps::kVerticalColumnCheckpointLayout)});
      }
      write_vertical_column_files(config, coordinate, result);
      const auto diagnostics = mps::diagnostics::diagnose_vertical_column(
          config, coordinate, result.state, result.mass_flux, result.maximum_cfl,
          result.budget);
      mps::write_run_metadata(std::cout, mps::make_run_metadata(config), config);
      write_vertical_column_result(std::cout, result, diagnostics);
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
      mps::ExperimentConfig run_config = config;
      std::optional<mps::ControlRequestV1> control_request;
      if (command_line.control_request_path.has_value()) {
        if (command_line.restart_path.has_value() ||
            command_line.checkpoint_path.has_value()) {
          throw std::invalid_argument(
              "machine control mode does not accept checkpoint or restart options");
        }
        control_request = mps::load_control_request(*command_line.control_request_path);
        run_config = mps::apply_control_request(config, *control_request);
        machine_events = std::make_unique<MachineEventWriter>();
        machine_events->run_id = control_request->run_id;
        machine_events->emit("run.accepted");
      }
      MachineFrameContext frame_context;
      mps::ShallowWaterRunHooks hooks;
      if (machine_events != nullptr) {
        machine_events->emit("run.started");
        frame_context = {.directory = run_config.output_directory,
                         .config_fingerprint = mps::config_fingerprint(run_config),
                         .cells_per_panel = run_config.grid.cells_per_panel,
                         .events = machine_events.get()};
        hooks = {.on_frame = write_machine_frame,
                 .observer_context = &frame_context,
                 .is_cancelled = cancellation_requested};
        std::signal(SIGINT, request_cancellation);
        std::signal(SIGTERM, request_cancellation);
      }
      const auto start = std::chrono::steady_clock::now();
      const auto result = mps::run_shallow_water(
          run_config, std::move(initial_state), command_line.stop_after_step,
          control_request.has_value() ? std::span<const mps::InitialConditionEditV1>(
                                            control_request->initial_edits)
                                      : std::span<const mps::InitialConditionEditV1>{},
          hooks);
      const mps::Real wall_time_s =
          std::chrono::duration<mps::Real>(std::chrono::steady_clock::now() - start)
              .count();
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
      write_shallow_water_snapshot(run_config, result);
      write_shallow_water_diagnostics(run_config, result);
      if (machine_events != nullptr) {
        std::ofstream metadata(
            std::filesystem::path(run_config.output_directory) / "run-metadata.txt",
            std::ios::trunc);
        if (!metadata) {
          throw std::runtime_error("unable to open machine run metadata");
        }
        mps::write_run_metadata(metadata, mps::make_run_metadata(run_config),
                                run_config);
        const auto& diagnostics = result.final_diagnostics;
        machine_events->emit("diagnostics.sample",
                             "\"timeSeconds\":" + std::to_string(result.state.time_s) +
                                 ",\"step\":" + std::to_string(result.state.step) +
                                 ",\"mass\":" + std::to_string(diagnostics.mass) +
                                 ",\"energy\":" + std::to_string(diagnostics.energy));
        if (g_cancel_requested.load() && !result.reached_end_time) {
          machine_events->emit("run.cancelled");
        } else if (result.reached_end_time) {
          machine_events->emit("run.completed");
        } else {
          machine_events->emit("run.cancelled");
        }
      } else {
        mps::write_run_metadata(std::cout, mps::make_run_metadata(run_config),
                                run_config);
        write_shallow_water_result(std::cout, result, wall_time_s, grid.cell_count());
        write_shallow_water_errors(std::cout, run_config, grid, result);
      }
    }
    return 0;
  } catch (const std::exception& error) {
    if (machine_events != nullptr) {
      try {
        machine_events->emit("run.failed", "\"code\":\"native_error\",\"message\":\"" +
                                               json_escape(error.what()) + "\"");
      } catch (...) {
      }
    }
    std::cerr << "error: " << error.what() << '\n';
    return 1;
  }
}
