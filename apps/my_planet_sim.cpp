#include <atomic>
#include <charconv>
#include <chrono>
#include <cmath>
#include <csignal>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <numbers>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "myplanetsim/config/experiment_config.hpp"
#include "myplanetsim/control/control_request.hpp"
#include "myplanetsim/diagnostics/climate_statistics.hpp"
#include "myplanetsim/diagnostics/dry_hydrostatic_diagnostics.hpp"
#include "myplanetsim/diagnostics/reductions.hpp"
#include "myplanetsim/diagnostics/vertical_column_diagnostics.hpp"
#include "myplanetsim/dynamics/dry_hydrostatic_driver.hpp"
#include "myplanetsim/dynamics/dry_hydrostatic_sources.hpp"
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

// Dry hydrostatic runs publish a FrameV2 snapshot instead of the shallow-water FrameV1
// payload, so the context carries the level count and the frame interval the control
// request asked for. The driver observes every accepted step; only the requested
// interval and the final step are published.
struct MachineFrameV2Context {
  std::filesystem::path directory;
  std::string config_fingerprint;
  mps::Index cells_per_panel = 0;
  mps::Index levels = 0;
  std::uint64_t frame_interval_steps = 1;
  MachineEventWriter* events = nullptr;
  std::uint64_t frame_sequence = 0;
  std::uint64_t last_published_step = std::numeric_limits<std::uint64_t>::max();
};

void write_machine_frame_v2(MachineFrameV2Context& frame,
                            const mps::DryHydrostaticState& state,
                            const mps::DryHydrostaticDerived& derived) {
  const auto sequence = frame.frame_sequence++;
  const std::string filename = "frame_" + std::to_string(sequence) + ".bin";
  const auto path = frame.directory / filename;
  mps::write_frame_v2_file(path,
                           {.cells_per_panel = frame.cells_per_panel,
                            .levels = frame.levels,
                            .time_s = state.time_s,
                            .step = state.step,
                            .config_fingerprint = frame.config_fingerprint},
                           state.surface_pressure_pa, derived);
  frame.last_published_step = state.step;
  frame.events->emit(
      "frame.ready",
      "\"frameSequence\":" + std::to_string(sequence) + ",\"timeSeconds\":" +
          std::to_string(state.time_s) + ",\"step\":" + std::to_string(state.step) +
          ",\"relativePath\":\"" + json_escape(filename) +
          "\",\"byteLength\":" + std::to_string(std::filesystem::file_size(path)) +
          ",\"frameSchemaVersion\":2");
}

void print_control_description() {
  std::cout << "{\"protocolVersion\":1,\"supportedEdits\":[\"gaussian_depth\"],"
               "\"maxEditCount\":64,\"maxEndTimeSeconds\":31536000,"
               "\"maxTimeStepSeconds\":86400,\"maxFrameIntervalSteps\":1000000,"
               "\"modelKinds\":[\"shallow_water\",\"dry_hydrostatic\"],"
               "\"frameSchemaVersions\":[1,2],"
               "\"dryHydrostatic\":{\"maxCellsPerPanel\":24,\"maxLevels\":30,"
               "\"supportedEdits\":[]}}\n";
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

struct PhysicsDiagnosticsRow {
  mps::Real time_s;
  std::uint64_t step;
  mps::HeldSuarezDiagnostics rates;
  mps::Real cumulative_thermal_energy_j;
  mps::Real cumulative_drag_energy_j;
  mps::Real measured_energy_change_j;
  mps::Real energy_residual_j;
  mps::DryHydrostaticDiagnostics health;
  mps::Real maximum_wind_m_s;
  std::uint64_t non_finite_count;
};

class PhysicsDiagnosticsAccumulator {
 public:
  PhysicsDiagnosticsAccumulator(const mps::ExperimentConfig& config,
                                const mps::CubedSphereGrid& grid)
      : config_(config), grid_(grid) {}

  void observe_step(const mps::DryHydrostaticState& state,
                    const mps::DryHydrostaticStepDiagnostics& step) {
    constexpr mps::Real spin_up_end_s = 200.0 * 86400.0;
    if (previous_time_s_.has_value() && *previous_time_s_ >= spin_up_end_s) {
      cumulative_thermal_energy_j_ += step.thermal_energy_contribution_j;
      cumulative_drag_energy_j_ += step.rayleigh_drag_energy_contribution_j;
    }
    previous_time_s_ = state.time_s;
  }

  void observe_sample(const mps::DryHydrostaticState& state,
                      const mps::DryHydrostaticDerived& derived,
                      const mps::DryHydrostaticStepDiagnostics& step) {
    constexpr mps::Real spin_up_end_s = 200.0 * 86400.0;
    const auto health =
        mps::diagnose_dry_hydrostatic_budgets(grid_, state, derived, config_.planet);
    if (state.time_s >= spin_up_end_s && !window_energy_j_.has_value()) {
      window_energy_j_ = health.total_energy_j;
    }

    mps::Real maximum_wind = 0.0;
    std::uint64_t non_finite_count = 0;
    for (const auto value : derived.temperature_k)
      non_finite_count += !std::isfinite(value);
    for (const auto value : derived.velocity_m_s) {
      non_finite_count += !mps::is_finite(value);
      if (mps::is_finite(value))
        maximum_wind = std::max(maximum_wind, mps::norm(value));
    }
    for (const auto value : state.surface_pressure_pa)
      non_finite_count += !std::isfinite(value);
    const mps::Real measured_energy_change =
        window_energy_j_.has_value() ? health.total_energy_j - *window_energy_j_ : 0.0;
    const mps::Real physics_energy =
        cumulative_thermal_energy_j_ + cumulative_drag_energy_j_;
    last_ = {.time_s = state.time_s,
             .step = state.step,
             .rates = step.physics_rates,
             .cumulative_thermal_energy_j = cumulative_thermal_energy_j_,
             .cumulative_drag_energy_j = cumulative_drag_energy_j_,
             .measured_energy_change_j = measured_energy_change,
             .energy_residual_j = measured_energy_change - physics_energy,
             .health = health,
             .maximum_wind_m_s = maximum_wind,
             .non_finite_count = non_finite_count};
    if (state.step % config_.diagnostics.interval_steps == 0) rows_.push_back(*last_);
  }

  void write() {
    if (last_.has_value() && (rows_.empty() || rows_.back().step != last_->step)) {
      rows_.push_back(*last_);
    }
    const std::filesystem::path directory(config_.output_directory);
    std::filesystem::create_directories(directory);
    std::ofstream output(directory / "physics_diagnostics.csv", std::ios::trunc);
    if (!output) throw std::runtime_error("unable to open physics diagnostics CSV");
    output << "time_s,step,potential_temperature_mass_rate_k_kg_s,"
              "eastward_momentum_rate_n,northward_momentum_rate_n,"
              "thermal_energy_rate_w,rayleigh_drag_work_w,total_physics_energy_rate_w,"
              "cumulative_thermal_energy_j,cumulative_drag_energy_j,"
              "cumulative_physics_energy_j,measured_energy_change_j,energy_residual_j,"
              "dry_mass_kg,tracer_mass_kg,minimum_temperature_k,maximum_wind_m_s,"
              "non_finite_count\n";
    output << std::setprecision(std::numeric_limits<mps::Real>::max_digits10);
    for (const auto& row : rows_) {
      const auto physics_rate =
          row.rates.thermal_energy_rate_w + row.rates.rayleigh_drag_work_w;
      const auto physics_energy =
          row.cumulative_thermal_energy_j + row.cumulative_drag_energy_j;
      output << row.time_s << ',' << row.step << ','
             << row.rates.potential_temperature_mass_rate_k_kg_s << ','
             << row.rates.eastward_momentum_rate_n << ','
             << row.rates.northward_momentum_rate_n << ','
             << row.rates.thermal_energy_rate_w << ',' << row.rates.rayleigh_drag_work_w
             << ',' << physics_rate << ',' << row.cumulative_thermal_energy_j << ','
             << row.cumulative_drag_energy_j << ',' << physics_energy << ','
             << row.measured_energy_change_j << ',' << row.energy_residual_j << ','
             << row.health.dry_mass_kg << ',' << row.health.tracer_mass_kg << ','
             << row.health.minimum_temperature_k << ',' << row.maximum_wind_m_s << ','
             << row.non_finite_count << '\n';
    }
    if (!output)
      throw std::runtime_error("failed while writing physics diagnostics CSV");
  }

 private:
  const mps::ExperimentConfig& config_;
  const mps::CubedSphereGrid& grid_;
  std::optional<mps::Real> previous_time_s_;
  std::optional<mps::Real> window_energy_j_;
  mps::Real cumulative_thermal_energy_j_ = 0.0;
  mps::Real cumulative_drag_energy_j_ = 0.0;
  std::optional<PhysicsDiagnosticsRow> last_;
  std::vector<PhysicsDiagnosticsRow> rows_;
};

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
    } else if (config.kind == mps::ExperimentKind::kDryHydrostatic) {
      if (command_line.integrator != mps::IntegratorKind::kSspRk3) {
        throw std::invalid_argument("dry hydrostatic supports only ssprk3");
      }
      if (command_line.control_request_path.has_value()) {
        if (config.physics.kind == mps::PhysicsKind::kSurfaceEnergyBalance) {
          throw std::invalid_argument(
              "machine control mode does not support surface energy balance");
        }
        if (command_line.restart_path.has_value() ||
            command_line.checkpoint_path.has_value()) {
          throw std::invalid_argument(
              "machine control mode does not accept checkpoint or restart options");
        }
        const auto control_request =
            mps::load_control_request(*command_line.control_request_path);
        if (!control_request.initial_edits.empty()) {
          throw std::invalid_argument(
              "dry hydrostatic runs do not accept initial condition edits");
        }
        const auto run_config = mps::apply_control_request(config, control_request);
        machine_events = std::make_unique<MachineEventWriter>();
        machine_events->run_id = control_request.run_id;
        machine_events->emit("run.accepted");
        const mps::DryHydrostaticDriver driver(run_config);
        auto state = driver.initial_state();
        MachineFrameV2Context frame_context{
            .directory = run_config.output_directory,
            .config_fingerprint = mps::config_fingerprint(run_config),
            .cells_per_panel = run_config.grid.cells_per_panel,
            .levels = run_config.vertical.levels,
            .frame_interval_steps = run_config.diagnostics.interval_steps,
            .events = machine_events.get()};
        machine_events->emit("run.started");
        std::signal(SIGINT, request_cancellation);
        std::signal(SIGTERM, request_cancellation);
        driver.advance(
            state, run_config.run.end_time_s,
            [&frame_context](const mps::DryHydrostaticState& sampled,
                             const mps::DryHydrostaticDerived* derived,
                             const mps::DryHydrostaticStepDiagnostics&) {
              if (derived != nullptr) {
                write_machine_frame_v2(frame_context, sampled, *derived);
              }
            },
            [] { return g_cancel_requested.load(); });
        const auto derived = driver.diagnose(state);
        if (frame_context.last_published_step != state.step) {
          write_machine_frame_v2(frame_context, state, derived);
        }
        std::ofstream metadata(
            std::filesystem::path(run_config.output_directory) / "run-metadata.txt",
            std::ios::trunc);
        if (!metadata) {
          throw std::runtime_error("unable to open machine run metadata");
        }
        mps::write_run_metadata(metadata, mps::make_run_metadata(run_config),
                                run_config);
        if (run_config.orography.kind != mps::OrographyKind::kFlat) {
          const auto sources =
              mps::dry_hydrostatic_sources(driver.grid(), derived, run_config.planet);
          const auto terrain = mps::diagnose_terrain_budgets(
              driver.grid(), state, derived, sources,
              driver.orography().surface_geopotential_m2_s2(), run_config.planet);
          mps::write_terrain_diagnostics(metadata, terrain);
        }
        const auto diagnostics = mps::diagnose_dry_hydrostatic_budgets(
            driver.grid(), state, derived, run_config.planet);
        machine_events->emit(
            "diagnostics.sample",
            "\"timeSeconds\":" + std::to_string(state.time_s) +
                ",\"step\":" + std::to_string(state.step) +
                ",\"mass\":" + std::to_string(diagnostics.dry_mass_kg) +
                ",\"energy\":" + std::to_string(diagnostics.total_energy_j));
        machine_events->emit(state.time_s < run_config.run.end_time_s
                                 ? "run.cancelled"
                                 : "run.completed");
        return 0;
      }
      mps::DryHydrostaticDriver driver(config);
      auto state = driver.initial_state();
      if (command_line.restart_path.has_value()) {
        const auto cells = driver.grid().cell_count();
        const auto levels = static_cast<std::size_t>(config.vertical.levels);
        const bool has_surface =
            config.physics.kind == mps::PhysicsKind::kSurfaceEnergyBalance;
        const auto layout = has_surface ? mps::kDryHydrostaticSurfaceCheckpointLayout
                                        : mps::kDryHydrostaticCheckpointLayout;
        const auto state_size = cells + 5 * cells * levels + (has_surface ? cells : 0);
        auto checkpoint = mps::read_checkpoint_file(*command_line.restart_path,
                                                    fingerprint, layout, state_size);
        state = has_surface ? mps::unflatten_dry_hydrostatic_surface_state(
                                  checkpoint.time_s, checkpoint.step, checkpoint.state,
                                  cells, levels)
                            : mps::unflatten_dry_hydrostatic_state(
                                  checkpoint.time_s, checkpoint.step, checkpoint.state,
                                  cells, levels);
      }
      std::optional<PhysicsDiagnosticsAccumulator> physics_diagnostics;
      std::optional<mps::ClimateStatisticsAccumulator> climate_statistics;
      if (config.physics.kind == mps::PhysicsKind::kHeldSuarez ||
          config.physics.kind == mps::PhysicsKind::kPlanetaryNewtonian) {
        physics_diagnostics.emplace(config, driver.grid());
        climate_statistics.emplace(driver.grid(),
                                   static_cast<std::size_t>(config.vertical.levels));
      }
      std::vector<std::pair<mps::Real, mps::SurfaceEnergyDiagnostics>>
          surface_diagnostics;
      driver.advance(
          state, config.run.end_time_s,
          [&physics_diagnostics, &climate_statistics, &surface_diagnostics, &config](
              const mps::DryHydrostaticState& sampled,
              const mps::DryHydrostaticDerived* derived,
              const mps::DryHydrostaticStepDiagnostics& step) {
            if (physics_diagnostics.has_value()) {
              physics_diagnostics->observe_step(sampled, step);
              if (derived != nullptr)
                physics_diagnostics->observe_sample(sampled, *derived, step);
            }
            if (climate_statistics.has_value() && derived != nullptr) {
              climate_statistics->observe(sampled, *derived);
            }
            if (config.physics.kind == mps::PhysicsKind::kSurfaceEnergyBalance &&
                derived != nullptr)
              surface_diagnostics.emplace_back(sampled.time_s, step.surface_rates);
          });
      if (physics_diagnostics.has_value()) physics_diagnostics->write();
      if (climate_statistics.has_value()) {
        const std::filesystem::path directory(config.output_directory);
        std::filesystem::create_directories(directory);
        std::ofstream output(directory / "climate_statistics.csv", std::ios::trunc);
        if (!output) throw std::runtime_error("unable to open climate statistics CSV");
        mps::write_climate_statistics_csv(output, climate_statistics->rows());
      }
      if (config.physics.kind == mps::PhysicsKind::kSurfaceEnergyBalance) {
        const std::filesystem::path directory(config.output_directory);
        std::filesystem::create_directories(directory);
        std::ofstream surface_state(directory / "surface_state.csv", std::ios::trunc);
        if (!surface_state)
          throw std::runtime_error("unable to open surface state CSV");
        surface_state << std::setprecision(17)
                      << "panel,i,j,longitude_deg,latitude_deg,height_m,"
                         "land_fraction,heat_capacity_j_m2_k,"
                         "surface_temperature_k\n";
        const auto& boundary = *driver.surface_boundary();
        for (std::size_t cell = 0; cell < driver.grid().cell_count(); ++cell) {
          const auto& geometry = driver.grid().cells()[cell];
          const auto longitude = std::atan2(geometry.center.y, geometry.center.x) *
                                 180.0 / std::numbers::pi_v<mps::Real>;
          const auto latitude = std::asin(std::clamp(geometry.center.z, -1.0, 1.0)) *
                                180.0 / std::numbers::pi_v<mps::Real>;
          const auto fraction = boundary.land_fraction()[cell];
          surface_state << mps::panel_index(geometry.id.panel) << ',' << geometry.id.i
                        << ',' << geometry.id.j << ',' << longitude << ',' << latitude
                        << ','
                        << boundary.surface_geopotential_m2_s2()[cell] /
                               config.planet.gravity_m_s2
                        << ',' << fraction << ','
                        << mps::mixed_surface_heat_capacity(
                               fraction, config.surface->land_heat_capacity_j_m2_k,
                               config.surface->ocean_heat_capacity_j_m2_k)
                        << ',' << state.surface_temperature_k[cell] << '\n';
        }
        std::ofstream diagnostics_output(directory / "surface_diagnostics.csv",
                                         std::ios::trunc);
        if (!diagnostics_output)
          throw std::runtime_error("unable to open surface diagnostics CSV");
        diagnostics_output
            << std::setprecision(17)
            << "time_s,absorbed_stellar_power_w,internal_heat_power_w,"
               "outgoing_longwave_power_w,sensible_to_atmosphere_power_w,"
               "surface_storage_rate_w,surface_budget_residual_w\n";
        for (const auto& [time, rates] : surface_diagnostics)
          diagnostics_output << time << ',' << rates.absorbed_stellar_power_w << ','
                             << rates.internal_heat_power_w << ','
                             << rates.outgoing_longwave_power_w << ','
                             << rates.sensible_to_atmosphere_power_w << ','
                             << rates.surface_storage_rate_w << ','
                             << rates.surface_budget_residual_w << '\n';
      }
      if (command_line.checkpoint_path.has_value()) {
        const bool has_surface =
            config.physics.kind == mps::PhysicsKind::kSurfaceEnergyBalance;
        mps::write_checkpoint_file(
            *command_line.checkpoint_path,
            {.time_s = state.time_s,
             .step = state.step,
             .state =
                 has_surface
                     ? mps::flatten_dry_hydrostatic_surface_state(
                           state, static_cast<std::size_t>(config.vertical.levels))
                     : mps::flatten_dry_hydrostatic_state(
                           state, static_cast<std::size_t>(config.vertical.levels)),
             .config_fingerprint = fingerprint,
             .layout_id =
                 std::string(has_surface ? mps::kDryHydrostaticSurfaceCheckpointLayout
                                         : mps::kDryHydrostaticCheckpointLayout)});
      }
      const auto derived = driver.diagnose(state);
      const auto diagnostics = mps::diagnose_dry_hydrostatic_budgets(
          driver.grid(), state, derived, config.planet);
      mps::write_run_metadata(std::cout, mps::make_run_metadata(config), config);
      if (config.orography.kind != mps::OrographyKind::kFlat) {
        const auto sources =
            mps::dry_hydrostatic_sources(driver.grid(), derived, config.planet);
        const auto terrain = mps::diagnose_terrain_budgets(
            driver.grid(), state, derived, sources,
            driver.orography().surface_geopotential_m2_s2(), config.planet);
        mps::write_terrain_diagnostics(std::cout, terrain);
      }
      std::cout << "result.status = complete\nresult.time_s = " << state.time_s
                << "\nresult.step = " << state.step
                << "\ndiagnostics.dry_mass_kg = " << diagnostics.dry_mass_kg
                << "\ndiagnostics.total_energy_j = " << diagnostics.total_energy_j
                << '\n';
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
