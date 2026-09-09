#include <algorithm>
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
#include <numbers>
#include <optional>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "myplanetsim/config/experiment_config.hpp"
#include "myplanetsim/diagnostics/climate_statistics.hpp"
#include "myplanetsim/diagnostics/dry_hydrostatic_diagnostics.hpp"
#include "myplanetsim/diagnostics/moist_diagnostics.hpp"
#include "myplanetsim/diagnostics/period_mean.hpp"
#include "myplanetsim/diagnostics/reductions.hpp"
#include "myplanetsim/diagnostics/vertical_column_diagnostics.hpp"
#include "myplanetsim/dynamics/dry_hydrostatic_driver.hpp"
#include "myplanetsim/dynamics/dry_hydrostatic_sources.hpp"
#include "myplanetsim/dynamics/shallow_water_benchmarks.hpp"
#include "myplanetsim/dynamics/shallow_water_driver.hpp"
#include "myplanetsim/dynamics/vertical_column_driver.hpp"
#include "myplanetsim/grid/cubed_sphere_grid.hpp"
#include "myplanetsim/io/checkpoint.hpp"
#include "myplanetsim/io/run_metadata.hpp"
#include "myplanetsim/io/visual_dataset.hpp"
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
  double progress_interval_s = 10.0;
};

void print_usage(std::ostream& output) {
  output << "Usage: my_planet_sim --config PATH [options]\n"
         << "Options:\n"
         << "  --integrator euler|ssprk3  Time integrator (default: ssprk3)\n"
         << "  --checkpoint PATH          Write the final or stopped state\n"
         << "  --restart PATH             Continue from a checkpoint\n"
         << "  --stop-after-step N        Stop after absolute step N\n"
         << "  --progress-interval-s SEC  Dry-core progress interval; 0 disables "
            "(default: 10)\n"
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

[[nodiscard]] double parse_progress_interval(const std::string_view text) {
  double value = 0.0;
  const auto result = std::from_chars(text.data(), text.data() + text.size(), value);
  if (result.ec != std::errc{} || result.ptr != text.data() + text.size() ||
      !std::isfinite(value) || value < 0.0) {
    throw std::invalid_argument("invalid --progress-interval-s value");
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
    } else if (argument == "--progress-interval-s") {
      command_line.progress_interval_s = parse_progress_interval(value);
    } else {
      throw std::invalid_argument("unknown option: " + std::string(argument));
    }
  }
  if (command_line.config_path.empty()) {
    throw std::invalid_argument("--config is required");
  }
  return command_line;
}

std::atomic_bool g_cancel_requested = false;

void request_cancellation(const int) noexcept { g_cancel_requested.store(true); }

class ProgressReporter {
 public:
  ProgressReporter(const mps::ExperimentConfig& config,
                   const mps::DryHydrostaticState& initial_state,
                   const double interval_s)
      : start_time_s_(config.run.start_time_s),
        end_time_s_(config.run.end_time_s),
        initial_model_time_s_(initial_state.time_s),
        interval_s_(interval_s),
        wall_start_(std::chrono::steady_clock::now()),
        last_report_(wall_start_) {}

  void observe(const mps::DryHydrostaticState& state) {
    if (!(interval_s_ > 0.0)) return;
    const auto now = std::chrono::steady_clock::now();
    if (reported_ &&
        std::chrono::duration<double>(now - last_report_).count() < interval_s_) {
      return;
    }
    emit(state, "running", now);
    reported_ = true;
    last_report_ = now;
  }

  void finish(const mps::DryHydrostaticState& state, const std::string_view status) {
    if (!(interval_s_ > 0.0)) return;
    emit(state, status, std::chrono::steady_clock::now());
  }

 private:
  [[nodiscard]] static std::string format_duration(const double seconds) {
    if (!std::isfinite(seconds) || seconds < 0.0) return "unknown";
    const auto total_seconds = static_cast<std::uint64_t>(std::round(seconds));
    const auto hours = total_seconds / 3600;
    const auto minutes = (total_seconds % 3600) / 60;
    const auto remainder = total_seconds % 60;
    std::ostringstream result;
    result << hours << ':' << std::setfill('0') << std::setw(2) << minutes << ':'
           << std::setw(2) << remainder;
    return result.str();
  }

  void emit(const mps::DryHydrostaticState& state, const std::string_view status,
            const std::chrono::steady_clock::time_point now) const {
    constexpr double seconds_per_day = 86400.0;
    const double duration_s = end_time_s_ - start_time_s_;
    const double completed_s = state.time_s - start_time_s_;
    const double percent = duration_s > 0.0
                               ? 100.0 * std::clamp(completed_s / duration_s, 0.0, 1.0)
                               : 100.0;
    const double elapsed_s = std::chrono::duration<double>(now - wall_start_).count();
    const double advanced_s = state.time_s - initial_model_time_s_;
    const double remaining_s = std::max(0.0, end_time_s_ - state.time_s);

    std::ostringstream line;
    line << std::fixed << std::setprecision(2) << "progress status=" << status
         << " step=" << state.step
         << " model_day=" << (state.time_s - start_time_s_) / seconds_per_day << '/'
         << duration_s / seconds_per_day << " percent=" << percent
         << " elapsed=" << format_duration(elapsed_s) << " eta=";
    if (advanced_s > 0.0) {
      line << format_duration(elapsed_s * remaining_s / advanced_s);
    } else {
      line << "unknown";
    }
    std::cerr << line.str() << '\n' << std::flush;
  }

  double start_time_s_;
  double end_time_s_;
  double initial_model_time_s_;
  double interval_s_;
  std::chrono::steady_clock::time_point wall_start_;
  std::chrono::steady_clock::time_point last_report_;
  bool reported_ = false;
};

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
  mps::Real cumulative_diffusion_energy_j;
  mps::Real measured_energy_change_j;
  mps::Real energy_residual_j;
  mps::DryHydrostaticDiagnostics health;
  mps::Real maximum_wind_m_s;
  std::uint64_t non_finite_count;
};

struct SurfaceDiagnosticsRow {
  mps::Real time_s;
  mps::SurfaceEnergyDiagnostics rates;
  mps::SurfaceEnergyBudget interval_budget;
};

struct RadiationDiagnosticsRow {
  mps::Real time_s;
  std::uint64_t step;
  mps::Real segment_start_time_s;
  std::uint64_t segment_start_step;
  mps::GrayRadiationDiagnostics rates;
  mps::GrayRadiationBudget interval_budget;
  mps::Real requested_time_step_s;
  mps::Real accepted_time_step_s;
  mps::Real radiation_stable_time_step_s;
  std::size_t cfl_retry_count;
  std::size_t invariant_retry_count;
  std::size_t solver_retry_count;
  std::size_t radiation_column_call_count;
  mps::Real radiation_wall_seconds;
};

struct SemiImplicitDiagnosticsRow {
  mps::Real time_s;
  std::uint64_t step;
  mps::DryHydrostaticStepDiagnostics diagnostics;
};

struct DryMixingDiagnosticsRow {
  mps::Real time_s;
  std::uint64_t step;
  mps::Real segment_start_time_s;
  std::uint64_t segment_start_step;
  mps::DryHydrostaticStepDiagnostics diagnostics;
};

struct TracerDiagnosticsRow {
  mps::Real time_s = 0.0;
  std::uint64_t step = 0;
  std::size_t tracer = 0;
  mps::TracerGlobalDiagnostics instantaneous;
  mps::Real interval_mass_change_kg = 0.0;
  mps::Real interval_physics_change_kg = 0.0;
  mps::Real interval_transport_diffusion_residual_kg = 0.0;
};

struct MoistureDiagnosticsRow {
  mps::Real time_s = 0.0;
  std::uint64_t step = 0;
  mps::Real interval_seconds = 0.0;
  mps::MoistureGlobalDiagnostics instantaneous;
  mps::Real dry_energy_j = 0.0;
  mps::Real atmospheric_latent_energy_j = 0.0;
  mps::Real surface_energy_j = 0.0;
  mps::Real dry_latent_surface_energy_j = 0.0;
  mps::MoistPhysicsStepDiagnostics interval;
  mps::Real cumulative_convective_precipitation_kg = 0.0;
  mps::Real cumulative_grid_scale_precipitation_kg = 0.0;
  mps::Real cumulative_evaporation_kg = 0.0;
  mps::Real cumulative_runoff_kg = 0.0;
  mps::Real cumulative_ocean_water_change_kg = 0.0;
  mps::Real cumulative_external_outflow_kg = 0.0;
  std::size_t physics_substep_count = 0;
  std::size_t physics_retry_count = 0;
};

void write_dry_mixing_diagnostics(const mps::ExperimentConfig& config,
                                  const mps::CubedSphereGrid& grid,
                                  const std::span<const DryMixingDiagnosticsRow> rows) {
  const std::filesystem::path directory(config.output_directory);
  std::filesystem::create_directories(directory);
  if (config.boundary_layer.kind != mps::BoundaryLayerKind::kNone) {
    std::ofstream output(directory / "boundary_layer_diagnostics.csv", std::ios::trunc);
    if (!output)
      throw std::runtime_error("unable to open boundary-layer diagnostics CSV");
    output << "time_s,step,segment_start_time_s,segment_start_step,requested_dt_s,"
              "accepted_dt_s,sensible_w_m2,surface_stress_impulse_n_s,mean_h_m,"
              "maximum_k_m_m2_s,maximum_k_h_m2_s,shallow_unresolved_area_fraction,"
              "model_top_area_fraction,atmospheric_heat_change_j,"
              "surface_heat_change_j,physical_shear_dissipation_j,"
              "physical_surface_drag_dissipation_j,backward_euler_dissipation_j,"
              "returned_dissipation_heat_j,heat_budget_residual_j,"
              "kinetic_energy_change_j,kinetic_energy_identity_residual_j,"
              "momentum_budget_residual_n_s,tracer_mass_change_kg,"
              "column_call_count,retry_count\n"
           << std::setprecision(std::numeric_limits<mps::Real>::max_digits10);
    for (const auto& row : rows) {
      const auto& step = row.diagnostics;
      const auto& value = step.boundary_layer;
      const mps::Real sensible_w_m2 =
          value.sensible_to_atmosphere_energy_j /
          (step.accepted_time_step_s * grid.total_area_m2());
      output << row.time_s << ',' << row.step << ',' << row.segment_start_time_s << ','
             << row.segment_start_step << ',' << step.requested_time_step_s << ','
             << step.accepted_time_step_s << ',' << sensible_w_m2 << ','
             << value.surface_stress_impulse_magnitude_n_s << ','
             << value.mean_boundary_layer_height_m << ','
             << value.maximum_momentum_diffusivity_m2_s << ','
             << value.maximum_heat_diffusivity_m2_s << ','
             << value.shallow_unresolved_area_fraction << ','
             << value.model_top_area_fraction << ',' << value.atmospheric_heat_change_j
             << ',' << value.surface_heat_change_j << ','
             << value.physical_shear_dissipation_j << ','
             << value.physical_surface_drag_dissipation_j << ','
             << value.backward_euler_dissipation_j << ','
             << value.returned_dissipation_heat_j << ',' << value.heat_budget_residual_j
             << ',' << value.kinetic_energy_change_j << ','
             << value.kinetic_energy_identity_residual_j << ','
             << value.momentum_budget_residual_n_s << ',' << value.tracer_mass_change_kg
             << ',' << step.boundary_layer_column_call_count << ',' << step.retry_count
             << '\n';
    }
  }
  if (config.convection.kind != mps::ConvectionKind::kNone) {
    std::ofstream output(directory / "convection_diagnostics.csv", std::ios::trunc);
    if (!output) throw std::runtime_error("unable to open convection diagnostics CSV");
    output << "time_s,step,segment_start_time_s,segment_start_step,requested_dt_s,"
              "accepted_dt_s,minimum_theta_difference_before_k,"
              "minimum_theta_difference_after_k,unstable_interface_fraction_before,"
              "unstable_interface_fraction_after,adjusted_column_count,"
              "adjusted_layer_count,adjusted_block_count,"
              "maximum_temperature_increment_k,enthalpy_change_j,"
              "dry_energy_attributed_change_j,column_call_count,retry_count\n"
           << std::setprecision(std::numeric_limits<mps::Real>::max_digits10);
    for (const auto& row : rows) {
      const auto& step = row.diagnostics;
      const auto& value = step.convection;
      output << row.time_s << ',' << row.step << ',' << row.segment_start_time_s << ','
             << row.segment_start_step << ',' << step.requested_time_step_s << ','
             << step.accepted_time_step_s << ','
             << value.minimum_theta_difference_before_k << ','
             << value.minimum_theta_difference_after_k << ','
             << value.unstable_interface_fraction_before << ','
             << value.unstable_interface_fraction_after << ','
             << value.adjusted_column_count << ',' << value.adjusted_layer_count << ','
             << value.adjusted_block_count << ','
             << value.maximum_temperature_increment_k << ',' << value.enthalpy_change_j
             << ',' << value.dry_energy_attributed_change_j << ','
             << step.convection_column_call_count << ',' << step.retry_count << '\n';
    }
  }
}

void write_semi_implicit_diagnostics(
    const mps::ExperimentConfig& config,
    const std::span<const SemiImplicitDiagnosticsRow> rows) {
  const std::filesystem::path directory(config.output_directory);
  std::filesystem::create_directories(directory);
  std::ofstream output(directory / "semi_implicit_diagnostics.csv", std::ios::trunc);
  if (!output) throw std::runtime_error("unable to open semi-implicit diagnostics CSV");
  output << "time_s,step,requested_dt_s,accepted_dt_s,advective_cfl,"
            "implicit_wave_cfl,vertical_cfl,selected_implicit_modes,"
            "linear_iterations_total,linear_iterations_maximum,"
            "linear_relative_residual_maximum,nonlinear_iterations,"
            "nonlinear_relative_residual,retry_count,cfl_retry_count,"
            "invariant_retry_count,solver_retry_count,wall_seconds_rhs,"
            "wall_seconds_linear_solve,wall_seconds_total\n"
         << std::setprecision(std::numeric_limits<mps::Real>::max_digits10);
  for (const auto& row : rows) {
    const auto& value = row.diagnostics;
    output << row.time_s << ',' << row.step << ',' << value.requested_time_step_s << ','
           << value.accepted_time_step_s << ',' << value.advective_cfl << ','
           << value.implicit_wave_courant << ',' << value.vertical_cfl << ','
           << value.selected_implicit_modes << ',' << value.linear_iterations_total
           << ',' << value.linear_iterations_maximum << ','
           << value.linear_relative_residual_maximum << ','
           << value.nonlinear_iterations << ',' << value.nonlinear_relative_residual
           << ',' << value.retry_count << ',' << value.cfl_retry_count << ','
           << value.invariant_retry_count << ',' << value.solver_retry_count << ','
           << value.wall_seconds_rhs << ',' << value.wall_seconds_linear_solve << ','
           << value.wall_seconds_total << '\n';
  }
  if (!output)
    throw std::runtime_error("failed while writing semi-implicit diagnostics CSV");
}

void add_surface_budget(mps::SurfaceEnergyBudget& total,
                        const mps::SurfaceEnergyBudget& step) {
  total.absorbed_stellar_energy_j += step.absorbed_stellar_energy_j;
  total.internal_heat_energy_j += step.internal_heat_energy_j;
  total.outgoing_longwave_energy_j += step.outgoing_longwave_energy_j;
  total.sensible_to_atmosphere_energy_j += step.sensible_to_atmosphere_energy_j;
  total.surface_storage_change_j += step.surface_storage_change_j;
  total.surface_budget_residual_j += step.surface_budget_residual_j;
}

void add_radiation_budget(mps::GrayRadiationBudget& total,
                          const mps::GrayRadiationBudget& step) {
  total.toa_incoming_shortwave_energy_j += step.toa_incoming_shortwave_energy_j;
  total.toa_reflected_shortwave_energy_j += step.toa_reflected_shortwave_energy_j;
  total.toa_outgoing_longwave_energy_j += step.toa_outgoing_longwave_energy_j;
  total.toa_net_upward_energy_j += step.toa_net_upward_energy_j;
  total.atmospheric_shortwave_heating_energy_j +=
      step.atmospheric_shortwave_heating_energy_j;
  total.atmospheric_longwave_heating_energy_j +=
      step.atmospheric_longwave_heating_energy_j;
  total.sensible_to_atmosphere_energy_j += step.sensible_to_atmosphere_energy_j;
  total.internal_heat_energy_j += step.internal_heat_energy_j;
  total.surface_storage_change_j += step.surface_storage_change_j;
  total.surface_time_integration_residual_j += step.surface_time_integration_residual_j;
  total.interface_conservation_residual_j += step.interface_conservation_residual_j;
  total.dry_thermal_energy_j += step.dry_thermal_energy_j;
  total.rayleigh_drag_energy_j += step.rayleigh_drag_energy_j;
}

void write_radiation_diagnostics(const mps::ExperimentConfig& config,
                                 const mps::CubedSphereGrid& grid,
                                 const std::span<const RadiationDiagnosticsRow> rows) {
  const std::filesystem::path directory(config.output_directory);
  std::filesystem::create_directories(directory);
  std::ofstream output(directory / "radiation_diagnostics.csv", std::ios::trunc);
  if (!output) throw std::runtime_error("unable to open radiation diagnostics CSV");
  output << "time_s,step,segment_start_time_s,segment_start_step,"
            "toa_incoming_sw_w_m2,toa_reflected_sw_w_m2,olr_w_m2,"
            "toa_net_upward_w_m2,surface_down_sw_w_m2,surface_up_sw_w_m2,"
            "surface_down_lw_w_m2,surface_up_lw_w_m2,atmospheric_sw_heating_w_m2,"
            "atmospheric_lw_heating_w_m2,surface_storage_w_m2,sensible_w_m2,"
            "internal_heat_w_m2,interface_conservation_residual_w_m2,"
            "dry_thermal_energy_rate_w,rayleigh_drag_work_w,"
            "interval_toa_incoming_sw_j,interval_toa_reflected_sw_j,"
            "interval_olr_j,interval_toa_net_upward_j,interval_atmospheric_sw_j,"
            "interval_atmospheric_lw_j,interval_sensible_j,interval_internal_heat_j,"
            "interval_surface_storage_change_j,"
            "interval_surface_time_integration_residual_j,"
            "interval_interface_conservation_residual_j,interval_dry_thermal_j,"
            "interval_drag_j,requested_dt_s,accepted_dt_s,"
            "radiation_stable_dt_s,cfl_retry_count,"
            "invariant_retry_count,solver_retry_count,radiation_column_call_count,"
            "radiation_wall_seconds\n"
         << std::setprecision(std::numeric_limits<mps::Real>::max_digits10);
  const mps::Real inverse_area = 1.0 / grid.total_area_m2();
  for (const auto& row : rows) {
    const auto& rate = row.rates;
    const auto& budget = row.interval_budget;
    output << row.time_s << ',' << row.step << ',' << row.segment_start_time_s << ','
           << row.segment_start_step << ','
           << rate.toa_incoming_shortwave_power_w * inverse_area << ','
           << rate.toa_reflected_shortwave_power_w * inverse_area << ','
           << rate.toa_outgoing_longwave_power_w * inverse_area << ','
           << rate.toa_net_upward_power_w * inverse_area << ','
           << rate.surface_down_shortwave_power_w * inverse_area << ','
           << rate.surface_up_shortwave_power_w * inverse_area << ','
           << rate.surface_down_longwave_power_w * inverse_area << ','
           << rate.surface_up_longwave_power_w * inverse_area << ','
           << rate.atmospheric_shortwave_heating_power_w * inverse_area << ','
           << rate.atmospheric_longwave_heating_power_w * inverse_area << ','
           << rate.surface_storage_rate_w * inverse_area << ','
           << rate.sensible_to_atmosphere_power_w * inverse_area << ','
           << rate.internal_heat_power_w * inverse_area << ','
           << rate.interface_conservation_residual_w * inverse_area << ','
           << rate.dry_thermal_energy_rate_w << ',' << rate.rayleigh_drag_work_w << ','
           << budget.toa_incoming_shortwave_energy_j << ','
           << budget.toa_reflected_shortwave_energy_j << ','
           << budget.toa_outgoing_longwave_energy_j << ','
           << budget.toa_net_upward_energy_j << ','
           << budget.atmospheric_shortwave_heating_energy_j << ','
           << budget.atmospheric_longwave_heating_energy_j << ','
           << budget.sensible_to_atmosphere_energy_j << ','
           << budget.internal_heat_energy_j << ',' << budget.surface_storage_change_j
           << ',' << budget.surface_time_integration_residual_j << ','
           << budget.interface_conservation_residual_j << ','
           << budget.dry_thermal_energy_j << ',' << budget.rayleigh_drag_energy_j << ','
           << row.requested_time_step_s << ',' << row.accepted_time_step_s << ','
           << row.radiation_stable_time_step_s << ',' << row.cfl_retry_count << ','
           << row.invariant_retry_count << ',' << row.solver_retry_count << ','
           << row.radiation_column_call_count << ',' << row.radiation_wall_seconds
           << '\n';
  }
  if (!output)
    throw std::runtime_error("failed while writing radiation diagnostics CSV");
}

void write_radiation_column(const mps::ExperimentConfig& config,
                            const mps::DryHydrostaticDriver& driver,
                            const mps::DryHydrostaticState& state,
                            const mps::DryHydrostaticDerived& derived) {
  const mps::AtmosphericHybridCoordinate coordinate(
      {config.vertical.a_half_pa, config.vertical.b_half},
      config.vertical.minimum_surface_pressure_pa,
      config.vertical.maximum_surface_pressure_pa,
      config.vertical.minimum_pressure_thickness_pa);
  const auto geometry = coordinate.geometry(
      state.surface_pressure_pa.front(), config.planet.gravity_m_s2,
      config.planet.gas_constant_j_kg_k, config.planet.heat_capacity_cp_j_kg_k,
      config.planet.reference_pressure_pa);
  const std::size_t levels = coordinate.levels();
  const auto orbit =
      mps::evaluate_orbit(*config.orbit, config.planet.rotation_rate_rad_s,
                          state.time_s - config.run.start_time_s);
  const mps::GrayRadiationColumnInput input{
      .pressure_half_pa = geometry.pressure_half_pa,
      .temperature_k = std::span<const mps::Real>(derived.temperature_k.data(), levels),
      .surface_temperature_k = state.surface_temperature_k.front(),
      .gravity_m_s2 = config.planet.gravity_m_s2,
      .stellar_flux_w_m2 = orbit.stellar_flux_w_m2,
      .cosine_solar_zenith =
          mps::cosine_solar_zenith(driver.grid().cells().front().center, orbit),
      .surface_albedo = config.surface->albedo,
      .surface_emissivity = config.surface->emissivity,
      .parameters = *config.radiation,
  };
  const auto column = mps::gray_radiation_column(input);
  const std::filesystem::path directory(config.output_directory);
  std::filesystem::create_directories(directory);
  std::ofstream output(directory / "radiation_column.csv", std::ios::trunc);
  if (!output) throw std::runtime_error("unable to open radiation column CSV");
  output << "time_s,cell_id,record_kind,index,pressure_pa,shortwave_optical_depth,"
            "longwave_optical_depth,shortwave_down_w_m2,shortwave_up_w_m2,"
            "longwave_down_w_m2,longwave_up_w_m2,net_flux_w_m2,temperature_k,"
            "shortwave_heating_w_m2,longwave_heating_w_m2,total_heating_w_m2\n"
         << std::setprecision(std::numeric_limits<mps::Real>::max_digits10);
  for (std::size_t interface = 0; interface <= levels; ++interface)
    output << state.time_s << ",0,interface," << interface << ','
           << geometry.pressure_half_pa[interface] << ",,,"
           << column.shortwave_down_w_m2[interface] << ','
           << column.shortwave_up_w_m2[interface] << ','
           << column.longwave_down_w_m2[interface] << ','
           << column.longwave_up_w_m2[interface] << ','
           << column.net_flux_w_m2[interface] << ",,,,\n";
  for (std::size_t level = 0; level < levels; ++level)
    output << state.time_s << ",0,full," << level << ','
           << geometry.pressure_full_pa[level] << ','
           << column.optical_depth.shortwave[level] << ','
           << column.optical_depth.longwave[level] << ",,,,,,"
           << derived.temperature_k[level] << ','
           << column.shortwave_convergence_w_m2[level] << ','
           << column.longwave_convergence_w_m2[level] << ','
           << column.radiative_convergence_w_m2[level] << '\n';
  if (!output) throw std::runtime_error("failed while writing radiation column CSV");
}

void write_mixing_column(const mps::ExperimentConfig& config,
                         const mps::DryHydrostaticDriver& driver,
                         const mps::DryHydrostaticState& state,
                         const mps::DryHydrostaticDerived& derived) {
  const std::size_t levels = static_cast<std::size_t>(config.vertical.levels);
  const mps::AtmosphericHybridCoordinate coordinate(
      {config.vertical.a_half_pa, config.vertical.b_half},
      config.vertical.minimum_surface_pressure_pa,
      config.vertical.maximum_surface_pressure_pa,
      config.vertical.minimum_pressure_thickness_pa);
  const auto geometry = coordinate.geometry(
      state.surface_pressure_pa.front(), config.planet.gravity_m_s2,
      config.planet.gas_constant_j_kg_k, config.planet.heat_capacity_cp_j_kg_k,
      config.planet.reference_pressure_pa);
  const auto theta =
      std::span<const mps::Real>(derived.potential_temperature_k.data(), levels);
  const auto temperature =
      std::span<const mps::Real>(derived.temperature_k.data(), levels);
  const auto velocity = std::span<const mps::Vec3>(derived.velocity_m_s.data(), levels);
  const auto tracer =
      std::span<const mps::Real>(derived.tracer_mixing_ratio.data(), levels);
  const auto& boundary = *driver.surface_boundary();
  const mps::Real surface_geopotential = boundary.surface_geopotential_m2_s2().front();
  const auto hydrostatic = mps::integrate_hydrostatic_column(
      geometry, std::vector<mps::Real>(theta.begin(), theta.end()),
      config.planet.heat_capacity_cp_j_kg_k, config.planet.gravity_m_s2,
      surface_geopotential);
  const mps::Real surface_height = surface_geopotential / config.planet.gravity_m_s2;
  std::vector<mps::Real> height_half(levels + 1);
  std::vector<mps::Real> height_full(levels);
  for (std::size_t interface = 0; interface <= levels; ++interface)
    height_half[interface] = hydrostatic.height_half_m[interface] - surface_height;
  height_half.back() = 0.0;
  for (std::size_t level = 0; level < levels; ++level)
    height_full[level] = hydrostatic.height_full_m[level] - surface_height;
  const auto bulk = mps::diagnose_boundary_layer_column(
      {.potential_temperature_k = theta,
       .temperature_k = temperature,
       .velocity_m_s = velocity,
       .pressure_half_pa = geometry.pressure_half_pa,
       .exner_half = geometry.exner_half,
       .height_half_m = height_half,
       .height_full_m = height_full,
       .surface_temperature_k = state.surface_temperature_k.front(),
       .surface_exner = geometry.exner_half.back(),
       .gravity_m_s2 = config.planet.gravity_m_s2,
       .gas_constant_j_kg_k = config.planet.gas_constant_j_kg_k,
       .heat_capacity_cp_j_kg_k = config.planet.heat_capacity_cp_j_kg_k,
       .critical_richardson = config.boundary_layer.critical_richardson,
       .turbulent_prandtl = config.boundary_layer.turbulent_prandtl,
       .gustiness_m_s = config.boundary_layer.gustiness_m_s,
       .land_fraction = boundary.land_fraction().front(),
       .land_roughness = {.momentum_m = config.surface->land_roughness_momentum_m,
                          .heat_m = config.surface->land_roughness_heat_m},
       .ocean_roughness = {.momentum_m = config.surface->ocean_roughness_momentum_m,
                           .heat_m = config.surface->ocean_roughness_heat_m}});

  const std::filesystem::path directory(config.output_directory);
  std::filesystem::create_directories(directory);
  std::ofstream output(directory / "mixing_column.csv", std::ios::trunc);
  if (!output) throw std::runtime_error("unable to open mixing column CSV");
  output << "time_s,cell_id,record_kind,index,pressure_pa,height_m,theta_k,"
            "temperature_k,velocity_x_m_s,velocity_y_m_s,velocity_z_m_s,tracer,"
            "density_kg_m3,k_m_m2_s,k_h_m2_s,k_q_m2_s,heat_flux_w_m2,"
            "momentum_flux_x_kg_m_s2,momentum_flux_y_kg_m_s2,"
            "momentum_flux_z_kg_m_s2,tracer_flux_kg_m2_s,boundary_layer_height_m,"
            "flux_kind\n"
         << std::setprecision(std::numeric_limits<mps::Real>::max_digits10);
  for (std::size_t interface = 0; interface <= levels; ++interface) {
    mps::Real heat_flux = 0.0;
    mps::Vec3 momentum_flux{};
    mps::Real tracer_flux = 0.0;
    std::string_view flux_kind = "diagnosed_zero_boundary";
    if (interface > 0 && interface < levels) {
      const mps::Real distance = height_full[interface - 1] - height_full[interface];
      heat_flux =
          bulk.density_half_kg_m3[interface] * config.planet.heat_capacity_cp_j_kg_k *
          geometry.exner_half[interface] * bulk.eddy_diffusivity_heat_m2_s[interface] *
          (theta[interface] - theta[interface - 1]) / distance;
      momentum_flux = bulk.density_half_kg_m3[interface] *
                      bulk.eddy_diffusivity_momentum_m2_s[interface] *
                      (velocity[interface] - velocity[interface - 1]) / distance;
      tracer_flux = bulk.density_half_kg_m3[interface] *
                    bulk.eddy_diffusivity_tracer_m2_s[interface] *
                    (tracer[interface] - tracer[interface - 1]) / distance;
      flux_kind = "final_diagnosed";
    } else if (interface == levels) {
      heat_flux = bulk.diagnostics.sensible_heat_flux_w_m2;
      momentum_flux = bulk.diagnostics.surface_stress_kg_m_s2;
      flux_kind = "final_diagnosed_surface";
    }
    output << state.time_s << ",0,interface," << interface << ','
           << geometry.pressure_half_pa[interface] << ',' << height_half[interface]
           << ",,,,,,," << bulk.density_half_kg_m3[interface] << ','
           << bulk.eddy_diffusivity_momentum_m2_s[interface] << ','
           << bulk.eddy_diffusivity_heat_m2_s[interface] << ','
           << bulk.eddy_diffusivity_tracer_m2_s[interface] << ',' << heat_flux << ','
           << momentum_flux.x << ',' << momentum_flux.y << ',' << momentum_flux.z << ','
           << tracer_flux << ',' << bulk.diagnostics.boundary_layer_height_m << ','
           << flux_kind << '\n';
  }
  for (std::size_t level = 0; level < levels; ++level)
    output << state.time_s << ",0,full," << level << ','
           << geometry.pressure_full_pa[level] << ',' << height_full[level] << ','
           << theta[level] << ',' << temperature[level] << ',' << velocity[level].x
           << ',' << velocity[level].y << ',' << velocity[level].z << ','
           << tracer[level] << ",,,,,,,,,," << bulk.diagnostics.boundary_layer_height_m
           << ",state_after_mixing\n";
  if (!output) throw std::runtime_error("failed while writing mixing column CSV");
}

void write_moist_column(const mps::ExperimentConfig& config,
                        const mps::DryHydrostaticDriver& driver,
                        const mps::DryHydrostaticState& state,
                        const mps::DryHydrostaticDerived& derived) {
  const std::size_t levels = static_cast<std::size_t>(config.vertical.levels);
  const std::size_t cells = driver.grid().cell_count();
  const mps::TracerRegistry registry(config.tracers);
  const std::size_t water = *registry.water_vapor_index();
  const mps::AtmosphericHybridCoordinate coordinate(
      {config.vertical.a_half_pa, config.vertical.b_half},
      config.vertical.minimum_surface_pressure_pa,
      config.vertical.maximum_surface_pressure_pa,
      config.vertical.minimum_pressure_thickness_pa);
  const auto geometry = coordinate.geometry(
      state.surface_pressure_pa.front(), config.planet.gravity_m_s2,
      config.planet.gas_constant_j_kg_k, config.planet.heat_capacity_cp_j_kg_k,
      config.planet.reference_pressure_pa);
  const auto theta =
      std::span<const mps::Real>(derived.potential_temperature_k.data(), levels);
  const auto hydrostatic = mps::integrate_hydrostatic_column(
      geometry, std::vector<mps::Real>(theta.begin(), theta.end()),
      config.planet.heat_capacity_cp_j_kg_k, config.planet.gravity_m_s2,
      driver.orography().surface_geopotential_m2_s2().front());
  const mps::DiluteMoistThermodynamics thermodynamics{
      .gas_constant_dry_air_j_kg_k = config.planet.gas_constant_j_kg_k,
      .heat_capacity_cp_j_kg_k = config.planet.heat_capacity_cp_j_kg_k};
  const std::filesystem::path directory(config.output_directory);
  std::filesystem::create_directories(directory);
  std::ofstream output(directory / "moist_column.csv", std::ios::trunc);
  if (!output) throw std::runtime_error("unable to open moist column CSV");
  output << "time_s,cell_id,level,pressure_pa,height_m,temperature_k,"
            "potential_temperature_k,vapor_mixing_ratio,saturation_mixing_ratio,"
            "relative_humidity,air_mass_kg_m2,surface_temperature_k,"
            "land_water_kg_m2,cumulative_convective_precipitation_kg,"
            "cumulative_grid_scale_precipitation_kg,cumulative_evaporation_kg,"
            "cumulative_runoff_kg\n"
         << std::setprecision(std::numeric_limits<mps::Real>::max_digits10);
  const mps::Real surface_height =
      driver.orography().surface_geopotential_m2_s2().front() /
      config.planet.gravity_m_s2;
  for (std::size_t level = 0; level < levels; ++level) {
    const std::size_t scalar = mps::dry_hydrostatic_offset(0, level, levels);
    const std::size_t vapor =
        mps::dry_hydrostatic_tracer_offset(water, 0, level, cells, levels);
    const mps::Real temperature = derived.temperature_k[scalar];
    const mps::Real q = derived.tracer_mixing_ratio[vapor];
    const mps::Real saturation = mps::saturation_mixing_ratio(
        temperature, derived.pressure_pa[scalar], thermodynamics);
    output << state.time_s << ",0," << level << ',' << derived.pressure_pa[scalar]
           << ',' << hydrostatic.height_full_m[level] - surface_height << ','
           << temperature << ',' << derived.potential_temperature_k[scalar] << ',' << q
           << ',' << saturation << ',' << q / saturation << ','
           << derived.air_mass_kg_m2[scalar] << ','
           << state.surface_temperature_k.front() << ','
           << state.land_water_kg_m2.front() << ','
           << state.cumulative_convective_precipitation_kg << ','
           << state.cumulative_grid_scale_precipitation_kg << ','
           << state.cumulative_evaporation_kg << ',' << state.cumulative_runoff_kg
           << '\n';
  }
  if (!output) throw std::runtime_error("failed while writing moist column CSV");

  std::ofstream surface_output(directory / "moist_surface_state.csv", std::ios::trunc);
  if (!surface_output)
    throw std::runtime_error("unable to open moist surface state CSV");
  surface_output << "time_s,cell_id,area_m2,land_fraction,surface_temperature_k,"
                    "land_water_kg_m2\n"
                 << std::setprecision(std::numeric_limits<mps::Real>::max_digits10);
  const auto& land_fraction = driver.surface_boundary()->land_fraction();
  for (std::size_t cell = 0; cell < cells; ++cell)
    surface_output << state.time_s << ',' << cell << ','
                   << driver.grid().cells()[cell].area_m2 << ',' << land_fraction[cell]
                   << ',' << state.surface_temperature_k[cell] << ','
                   << state.land_water_kg_m2[cell] << '\n';
  if (!surface_output)
    throw std::runtime_error("failed while writing moist surface state CSV");
}

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
      cumulative_diffusion_energy_j_ += step.diffusion_energy_contribution_j;
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
    const mps::Real attributed_energy = physics_energy + cumulative_diffusion_energy_j_;
    last_ = {.time_s = state.time_s,
             .step = state.step,
             .rates = step.physics_rates,
             .cumulative_thermal_energy_j = cumulative_thermal_energy_j_,
             .cumulative_drag_energy_j = cumulative_drag_energy_j_,
             .cumulative_diffusion_energy_j = cumulative_diffusion_energy_j_,
             .measured_energy_change_j = measured_energy_change,
             .energy_residual_j = measured_energy_change - attributed_energy,
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
              "cumulative_physics_energy_j,cumulative_diffusion_energy_j,"
              "cumulative_attributed_energy_j,measured_energy_change_j,"
              "energy_residual_j,"
              "dry_mass_kg,tracer_mass_kg,minimum_temperature_k,maximum_wind_m_s,"
              "non_finite_count\n";
    output << std::setprecision(std::numeric_limits<mps::Real>::max_digits10);
    for (const auto& row : rows_) {
      const auto physics_rate =
          row.rates.thermal_energy_rate_w + row.rates.rayleigh_drag_work_w;
      const auto physics_energy =
          row.cumulative_thermal_energy_j + row.cumulative_drag_energy_j;
      const auto attributed_energy = physics_energy + row.cumulative_diffusion_energy_j;
      output << row.time_s << ',' << row.step << ','
             << row.rates.potential_temperature_mass_rate_k_kg_s << ','
             << row.rates.eastward_momentum_rate_n << ','
             << row.rates.northward_momentum_rate_n << ','
             << row.rates.thermal_energy_rate_w << ',' << row.rates.rayleigh_drag_work_w
             << ',' << physics_rate << ',' << row.cumulative_thermal_energy_j << ','
             << row.cumulative_drag_energy_j << ',' << physics_energy << ','
             << row.cumulative_diffusion_energy_j << ',' << attributed_energy << ','
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
  mps::Real cumulative_diffusion_energy_j_ = 0.0;
  std::optional<PhysicsDiagnosticsRow> last_;
  std::vector<PhysicsDiagnosticsRow> rows_;
};

class MoistDiagnosticsAccumulator {
 public:
  MoistDiagnosticsAccumulator(const mps::ExperimentConfig& config,
                              const mps::CubedSphereGrid& grid,
                              const std::span<const mps::Real> land_fraction)
      : config_(config),
        grid_(grid),
        registry_(config.tracers),
        water_vapor_tracer_(*registry_.water_vapor_index()),
        thermodynamics_{
            .gas_constant_dry_air_j_kg_k = config.planet.gas_constant_j_kg_k,
            .heat_capacity_cp_j_kg_k = config.planet.heat_capacity_cp_j_kg_k},
        land_fraction_(land_fraction.begin(), land_fraction.end()) {
    if (land_fraction_.size() != grid_.cell_count())
      throw std::invalid_argument("moist diagnostics land-fraction shape mismatch");
  }

  void observe_step(const mps::DryHydrostaticStepDiagnostics& step) {
    if (!(step.accepted_time_step_s > 0.0)) return;
    pending_seconds_ += step.accepted_time_step_s;
    pending_substeps_ += step.physics_substep_count;
    pending_retries_ += step.physics_retry_count;
    const auto& local = step.moisture;
    pending_.evaporation_kg += local.evaporation_kg;
    pending_.convective_precipitation_kg += local.convective_precipitation_kg;
    pending_.grid_scale_precipitation_kg += local.grid_scale_precipitation_kg;
    pending_.runoff_kg += local.runoff_kg;
    pending_.ocean_water_change_kg += local.ocean_water_change_kg;
    pending_.external_outflow_kg += local.external_outflow_kg;
    pending_.water_budget_residual_kg += local.water_budget_residual_kg;
    pending_.moist_enthalpy_budget_residual_j += local.moist_enthalpy_budget_residual_j;
    pending_.maximum_relative_humidity =
        std::max(pending_.maximum_relative_humidity, local.maximum_relative_humidity);
    pending_.maximum_temperature_increment_k =
        std::max(pending_.maximum_temperature_increment_k,
                 local.maximum_temperature_increment_k);
    pending_.maximum_vapor_increment =
        std::max(pending_.maximum_vapor_increment, local.maximum_vapor_increment);
    pending_.cape_area_time_integral_j_m2_s_kg +=
        local.cape_area_time_integral_j_m2_s_kg;
    pending_.cin_area_time_integral_j_m2_s_kg += local.cin_area_time_integral_j_m2_s_kg;
    pending_.lcl_pressure_area_time_integral_pa_m2_s +=
        local.lcl_pressure_area_time_integral_pa_m2_s;
    pending_.convection_top_pressure_area_time_integral_pa_m2_s +=
        local.convection_top_pressure_area_time_integral_pa_m2_s;
    pending_.active_area_time_m2_s += local.active_area_time_m2_s;
    pending_.deep_area_time_m2_s += local.deep_area_time_m2_s;
    pending_.shallow_area_time_m2_s += local.shallow_area_time_m2_s;
    pending_.inactive_area_time_m2_s += local.inactive_area_time_m2_s;
    pending_.model_top_area_time_m2_s += local.model_top_area_time_m2_s;
    pending_.deep_column_count += local.deep_column_count;
    pending_.shallow_column_count += local.shallow_column_count;
    pending_.inactive_column_count += local.inactive_column_count;
  }

  void observe_sample(const mps::DryHydrostaticState& state,
                      const mps::DryHydrostaticDerived& derived) {
    const auto tracers = mps::diagnose_tracers(grid_, state, derived, registry_);
    const auto instantaneous = mps::diagnose_moisture(
        grid_, state, derived, land_fraction_, water_vapor_tracer_, thermodynamics_);
    const auto dry =
        mps::diagnose_dry_hydrostatic_budgets(grid_, state, derived, config_.planet);
    mps::Real surface_energy_j = 0.0;
    for (std::size_t cell = 0; cell < grid_.cell_count(); ++cell) {
      const mps::Real heat_capacity = mps::mixed_surface_heat_capacity(
          land_fraction_[cell], config_.surface->land_heat_capacity_j_m2_k,
          config_.surface->ocean_heat_capacity_j_m2_k);
      surface_energy_j += grid_.cells()[cell].area_m2 * heat_capacity *
                          state.surface_temperature_k[cell];
    }
    const mps::Real latent_energy_j = thermodynamics_.latent_heat_vaporization_j_kg *
                                      instantaneous.atmospheric_water_kg;
    if (previous_tracer_mass_kg_.empty()) {
      previous_tracer_mass_kg_.reserve(tracers.size());
      for (const auto& tracer : tracers)
        previous_tracer_mass_kg_.push_back(tracer.mass_kg);
    }
    const mps::Real water_physics_change = pending_.evaporation_kg -
                                           pending_.convective_precipitation_kg -
                                           pending_.grid_scale_precipitation_kg;
    for (std::size_t tracer = 0; tracer < tracers.size(); ++tracer) {
      const mps::Real mass_change =
          tracers[tracer].mass_kg - previous_tracer_mass_kg_[tracer];
      const mps::Real physics_change =
          tracer == water_vapor_tracer_ ? water_physics_change : 0.0;
      tracer_rows_.push_back(
          {.time_s = state.time_s,
           .step = state.step,
           .tracer = tracer,
           .instantaneous = tracers[tracer],
           .interval_mass_change_kg = mass_change,
           .interval_physics_change_kg = physics_change,
           .interval_transport_diffusion_residual_kg = mass_change - physics_change});
      previous_tracer_mass_kg_[tracer] = tracers[tracer].mass_kg;
    }
    moisture_rows_.push_back(
        {.time_s = state.time_s,
         .step = state.step,
         .interval_seconds = pending_seconds_,
         .instantaneous = instantaneous,
         .dry_energy_j = dry.total_energy_j,
         .atmospheric_latent_energy_j = latent_energy_j,
         .surface_energy_j = surface_energy_j,
         .dry_latent_surface_energy_j =
             dry.total_energy_j + latent_energy_j + surface_energy_j,
         .interval = pending_,
         .cumulative_convective_precipitation_kg =
             state.cumulative_convective_precipitation_kg,
         .cumulative_grid_scale_precipitation_kg =
             state.cumulative_grid_scale_precipitation_kg,
         .cumulative_evaporation_kg = state.cumulative_evaporation_kg,
         .cumulative_runoff_kg = state.cumulative_runoff_kg,
         .cumulative_ocean_water_change_kg = state.cumulative_ocean_water_change_kg,
         .cumulative_external_outflow_kg = state.cumulative_external_outflow_kg,
         .physics_substep_count = pending_substeps_,
         .physics_retry_count = pending_retries_});
    pending_ = {};
    pending_seconds_ = 0.0;
    pending_substeps_ = 0;
    pending_retries_ = 0;
  }

  void write() const {
    const std::filesystem::path directory(config_.output_directory);
    std::filesystem::create_directories(directory);
    std::ofstream tracer_output(directory / "tracer_diagnostics.csv", std::ios::trunc);
    if (!tracer_output)
      throw std::runtime_error("unable to open tracer diagnostics CSV");
    tracer_output
        << "time_s,step,tracer_index,tracer_name,tracer_role,mass_kg,"
           "minimum_mixing_ratio,maximum_mixing_ratio,interval_mass_change_kg,"
           "interval_physics_change_kg,"
           "interval_transport_diffusion_residual_kg\n"
        << std::setprecision(std::numeric_limits<mps::Real>::max_digits10);
    for (const auto& row : tracer_rows_) {
      const auto& descriptor = registry_[row.tracer];
      tracer_output << row.time_s << ',' << row.step << ',' << row.tracer << ','
                    << descriptor.name << ',' << mps::tracer_role_name(descriptor.role)
                    << ',' << row.instantaneous.mass_kg << ','
                    << row.instantaneous.minimum_mixing_ratio << ','
                    << row.instantaneous.maximum_mixing_ratio << ','
                    << row.interval_mass_change_kg << ','
                    << row.interval_physics_change_kg << ','
                    << row.interval_transport_diffusion_residual_kg << '\n';
    }

    std::ofstream moisture_output(directory / "moisture_diagnostics.csv",
                                  std::ios::trunc);
    if (!moisture_output)
      throw std::runtime_error("unable to open moisture diagnostics CSV");
    moisture_output
        << "time_s,step,interval_seconds,atmospheric_water_kg,"
           "precipitable_water_kg_m2,area_mean_relative_humidity,"
           "minimum_relative_humidity,maximum_relative_humidity,land_water_kg,"
           "area_mean_land_water_kg_m2,maximum_vapor_mixing_ratio,"
           "dilute_limit_exceedance_area_fraction,dry_energy_j,"
           "atmospheric_latent_energy_j,surface_energy_j,"
           "dry_latent_surface_energy_j,net_evaporation_kg,"
           "convective_precipitation_kg,grid_scale_precipitation_kg,runoff_kg,"
           "ocean_water_change_kg,external_outflow_kg,water_budget_residual_kg,"
           "moist_enthalpy_budget_residual_j,net_evaporation_kg_m2_s,"
           "precipitation_kg_m2_s,precipitation_mm_day_display,"
           "latent_heat_flux_w_m2,cumulative_evaporation_kg,"
           "cumulative_convective_precipitation_kg,"
           "cumulative_grid_scale_precipitation_kg,cumulative_runoff_kg,"
           "cumulative_ocean_water_change_kg,cumulative_external_outflow_kg,"
           "physics_substep_count,physics_retry_count\n"
        << std::setprecision(std::numeric_limits<mps::Real>::max_digits10);
    const mps::Real area = grid_.total_area_m2();
    for (const auto& row : moisture_rows_) {
      const mps::Real inverse_area_time =
          row.interval_seconds > 0.0 ? 1.0 / (area * row.interval_seconds) : 0.0;
      const mps::Real evaporation_rate =
          row.interval.evaporation_kg * inverse_area_time;
      const mps::Real precipitation_rate = (row.interval.convective_precipitation_kg +
                                            row.interval.grid_scale_precipitation_kg) *
                                           inverse_area_time;
      moisture_output
          << row.time_s << ',' << row.step << ',' << row.interval_seconds << ','
          << row.instantaneous.atmospheric_water_kg << ','
          << row.instantaneous.precipitable_water_kg_m2 << ','
          << row.instantaneous.area_mean_relative_humidity << ','
          << row.instantaneous.minimum_relative_humidity << ','
          << row.instantaneous.maximum_relative_humidity << ','
          << row.instantaneous.land_water_kg << ','
          << row.instantaneous.area_mean_land_water_kg_m2 << ','
          << row.instantaneous.maximum_vapor_mixing_ratio << ','
          << row.instantaneous.dilute_limit_exceedance_area_fraction << ','
          << row.dry_energy_j << ',' << row.atmospheric_latent_energy_j << ','
          << row.surface_energy_j << ',' << row.dry_latent_surface_energy_j << ','
          << row.interval.evaporation_kg << ','
          << row.interval.convective_precipitation_kg << ','
          << row.interval.grid_scale_precipitation_kg << ',' << row.interval.runoff_kg
          << ',' << row.interval.ocean_water_change_kg << ','
          << row.interval.external_outflow_kg << ','
          << row.interval.water_budget_residual_kg << ','
          << row.interval.moist_enthalpy_budget_residual_j << ',' << evaporation_rate
          << ',' << precipitation_rate << ',' << precipitation_rate * 86400.0 << ','
          << evaporation_rate * thermodynamics_.latent_heat_vaporization_j_kg << ','
          << row.cumulative_evaporation_kg << ','
          << row.cumulative_convective_precipitation_kg << ','
          << row.cumulative_grid_scale_precipitation_kg << ','
          << row.cumulative_runoff_kg << ',' << row.cumulative_ocean_water_change_kg
          << ',' << row.cumulative_external_outflow_kg << ','
          << row.physics_substep_count << ',' << row.physics_retry_count << '\n';
    }

    std::ofstream convection_output(directory / "moist_convection_diagnostics.csv",
                                    std::ios::trunc);
    if (!convection_output)
      throw std::runtime_error("unable to open moist convection diagnostics CSV");
    convection_output
        << "time_s,step,interval_seconds,area_mean_cape_j_kg,"
           "area_mean_cin_j_kg,active_mean_lcl_pressure_pa,"
           "active_mean_convection_top_pressure_pa,deep_area_fraction,"
           "shallow_area_fraction,inactive_area_fraction,model_top_area_fraction,"
           "maximum_temperature_increment_k,maximum_vapor_increment,"
           "moist_enthalpy_budget_residual_j,deep_column_calls,"
           "shallow_column_calls,inactive_column_calls,physics_substep_count,"
           "physics_retry_count\n"
        << std::setprecision(std::numeric_limits<mps::Real>::max_digits10);
    for (const auto& row : moisture_rows_) {
      const auto& value = row.interval;
      const mps::Real area_time = area * row.interval_seconds;
      const mps::Real inverse_area_time = area_time > 0.0 ? 1.0 / area_time : 0.0;
      const mps::Real inverse_active =
          value.active_area_time_m2_s > 0.0 ? 1.0 / value.active_area_time_m2_s : 0.0;
      const auto area_fraction = [inverse_area_time](const mps::Real value) {
        return std::clamp(value * inverse_area_time, 0.0, 1.0);
      };
      convection_output
          << row.time_s << ',' << row.step << ',' << row.interval_seconds << ','
          << value.cape_area_time_integral_j_m2_s_kg * inverse_area_time << ','
          << value.cin_area_time_integral_j_m2_s_kg * inverse_area_time << ','
          << value.lcl_pressure_area_time_integral_pa_m2_s * inverse_active << ','
          << value.convection_top_pressure_area_time_integral_pa_m2_s * inverse_active
          << ',' << area_fraction(value.deep_area_time_m2_s) << ','
          << area_fraction(value.shallow_area_time_m2_s) << ','
          << area_fraction(value.inactive_area_time_m2_s) << ','
          << area_fraction(value.model_top_area_time_m2_s) << ','
          << value.maximum_temperature_increment_k << ','
          << value.maximum_vapor_increment << ','
          << value.moist_enthalpy_budget_residual_j << ',' << value.deep_column_count
          << ',' << value.shallow_column_count << ',' << value.inactive_column_count
          << ',' << row.physics_substep_count << ',' << row.physics_retry_count << '\n';
    }
    if (!tracer_output || !moisture_output || !convection_output)
      throw std::runtime_error("failed while writing moist diagnostics CSV");
  }

 private:
  const mps::ExperimentConfig& config_;
  const mps::CubedSphereGrid& grid_;
  mps::TracerRegistry registry_;
  std::size_t water_vapor_tracer_ = 0;
  mps::DiluteMoistThermodynamics thermodynamics_;
  std::vector<mps::Real> land_fraction_;
  mps::MoistPhysicsStepDiagnostics pending_;
  mps::Real pending_seconds_ = 0.0;
  std::size_t pending_substeps_ = 0;
  std::size_t pending_retries_ = 0;
  std::vector<mps::Real> previous_tracer_mass_kg_;
  std::vector<TracerDiagnosticsRow> tracer_rows_;
  std::vector<MoistureDiagnosticsRow> moisture_rows_;
};

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
    } else if (config.kind == mps::ExperimentKind::kVerticalColumn) {
      if (command_line.integrator != mps::IntegratorKind::kSspRk3)
        throw std::invalid_argument("vertical column supports only ssprk3");
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
      mps::DryHydrostaticDriver driver(config);
      auto state = driver.initial_state();
      if (command_line.restart_path.has_value()) {
        const auto cells = driver.grid().cell_count();
        const auto levels = static_cast<std::size_t>(config.vertical.levels);
        const bool has_surface =
            config.physics.kind == mps::PhysicsKind::kSurfaceEnergyBalance ||
            config.physics.kind == mps::PhysicsKind::kGrayRadiation;
        const bool multitracer = !config.tracers.empty();
        const bool moist = config.moisture.kind == mps::MoistureKind::kDiluteWater;
        const mps::TracerRegistry registry = multitracer
                                                 ? mps::TracerRegistry(config.tracers)
                                                 : mps::TracerRegistry::legacy();
        const std::string layout =
            moist ? mps::dry_hydrostatic_moist_checkpoint_layout(registry)
            : multitracer
                ? mps::dry_hydrostatic_multitracer_checkpoint_layout(registry,
                                                                     has_surface)
                : std::string(has_surface ? mps::kDryHydrostaticSurfaceCheckpointLayout
                                          : mps::kDryHydrostaticCheckpointLayout);
        const auto state_size =
            moist ? cells + (4 + registry.size()) * cells * levels + 2 * cells + 6
            : multitracer ? cells + (4 + registry.size()) * cells * levels +
                                (has_surface ? cells : 0)
                          : cells + 5 * cells * levels + (has_surface ? cells : 0);
        auto checkpoint = mps::read_checkpoint_file(*command_line.restart_path,
                                                    fingerprint, layout, state_size);
        if (moist)
          state = mps::unflatten_dry_hydrostatic_moist_state(
              checkpoint.time_s, checkpoint.step, checkpoint.state, cells, levels,
              registry.size());
        else if (multitracer)
          state = mps::unflatten_dry_hydrostatic_multitracer_state(
              checkpoint.time_s, checkpoint.step, checkpoint.state, cells, levels,
              registry.size(), has_surface);
        else
          state = has_surface ? mps::unflatten_dry_hydrostatic_surface_state(
                                    checkpoint.time_s, checkpoint.step,
                                    checkpoint.state, cells, levels)
                              : mps::unflatten_dry_hydrostatic_state(
                                    checkpoint.time_s, checkpoint.step,
                                    checkpoint.state, cells, levels);
      }
      std::optional<PhysicsDiagnosticsAccumulator> physics_diagnostics;
      std::optional<mps::ClimateStatisticsAccumulator> climate_statistics;
      std::optional<MoistDiagnosticsAccumulator> moist_diagnostics;
      if (config.physics.kind == mps::PhysicsKind::kHeldSuarez ||
          config.physics.kind == mps::PhysicsKind::kPlanetaryNewtonian ||
          config.physics.kind == mps::PhysicsKind::kGrayRadiation) {
        physics_diagnostics.emplace(config, driver.grid());
        climate_statistics.emplace(driver.grid(),
                                   static_cast<std::size_t>(config.vertical.levels));
      }
      if (config.moisture.kind == mps::MoistureKind::kDiluteWater)
        moist_diagnostics.emplace(config, driver.grid(),
                                  driver.surface_boundary()->land_fraction());
      const auto visual_fields =
          mps::dry_visual_fields(config, driver.grid().cell_count(),
                                 static_cast<std::size_t>(config.vertical.levels),
                                 !state.surface_temperature_k.empty());
      const std::span<const mps::Real> visual_land_fraction =
          driver.surface_boundary().has_value()
              ? driver.surface_boundary()->land_fraction()
              : std::span<const mps::Real>{};
      mps::VisualDatasetWriter visual_writer(
          std::filesystem::path(config.output_directory) / "viewer", config,
          driver.grid(), driver.orography(), visual_land_fraction, visual_fields,
          fingerprint);
      std::optional<mps::diagnostics::PeriodMeanAccumulator> period_mean;
      if (config.statistics.enabled) {
        if (command_line.restart_path.has_value() &&
            std::filesystem::exists(
                mps::statistics_sidecar_path(*command_line.restart_path))) {
          const auto sidecar = mps::read_statistics_sidecar_file(
              mps::statistics_sidecar_path(*command_line.restart_path),
              mps::file_fnv1a64(*command_line.restart_path), state.time_s, state.step,
              visual_writer.value_count(), config.statistics.start_time_s,
              config.statistics.period_s);
          visual_writer.set_periods(sidecar.periods);
          period_mean.emplace(sidecar.accumulator);
        } else {
          period_mean.emplace(visual_writer.value_count(),
                              config.statistics.start_time_s,
                              config.statistics.period_s);
        }
      }
      std::vector<SurfaceDiagnosticsRow> surface_diagnostics;
      std::vector<RadiationDiagnosticsRow> radiation_diagnostics;
      std::vector<SemiImplicitDiagnosticsRow> semi_implicit_diagnostics;
      std::vector<DryMixingDiagnosticsRow> dry_mixing_diagnostics;
      mps::SurfaceEnergyBudget pending_surface_budget;
      mps::GrayRadiationBudget pending_radiation_budget;
      std::size_t pending_cfl_retries = 0;
      std::size_t pending_invariant_retries = 0;
      std::size_t pending_solver_retries = 0;
      std::size_t pending_radiation_column_calls = 0;
      mps::Real pending_radiation_wall_seconds = 0.0;
      const mps::Real segment_start_time_s = state.time_s;
      const std::uint64_t segment_start_step = state.step;
      ProgressReporter progress(config, state, command_line.progress_interval_s);
      g_cancel_requested.store(false);
      std::signal(SIGINT, request_cancellation);
      std::signal(SIGTERM, request_cancellation);
      driver.advance(
          state, config.run.end_time_s,
          [&physics_diagnostics, &climate_statistics, &moist_diagnostics,
           &surface_diagnostics, &radiation_diagnostics, &pending_surface_budget,
           &pending_radiation_budget, &pending_cfl_retries, &pending_invariant_retries,
           &pending_solver_retries, &pending_radiation_column_calls,
           &pending_radiation_wall_seconds, &semi_implicit_diagnostics,
           &dry_mixing_diagnostics, &progress, segment_start_time_s, segment_start_step,
           &config](const mps::DryHydrostaticState& sampled,
                    const mps::DryHydrostaticDerived* derived,
                    const mps::DryHydrostaticStepDiagnostics& step) {
            progress.observe(sampled);
            if (physics_diagnostics.has_value()) {
              physics_diagnostics->observe_step(sampled, step);
              if (derived != nullptr)
                physics_diagnostics->observe_sample(sampled, *derived, step);
            }
            if (climate_statistics.has_value() && derived != nullptr) {
              climate_statistics->observe(sampled, *derived);
            }
            if (moist_diagnostics.has_value()) {
              moist_diagnostics->observe_step(step);
              if (derived != nullptr)
                moist_diagnostics->observe_sample(sampled, *derived);
            }
            if (config.semi_implicit.has_value() && derived != nullptr) {
              semi_implicit_diagnostics.push_back({.time_s = sampled.time_s,
                                                   .step = sampled.step,
                                                   .diagnostics = step});
            }
            if ((config.convection.kind != mps::ConvectionKind::kNone ||
                 config.boundary_layer.kind != mps::BoundaryLayerKind::kNone) &&
                step.accepted_time_step_s > 0.0) {
              dry_mixing_diagnostics.push_back(
                  {.time_s = sampled.time_s,
                   .step = sampled.step,
                   .segment_start_time_s = segment_start_time_s,
                   .segment_start_step = segment_start_step,
                   .diagnostics = step});
            }
            if (config.physics.kind == mps::PhysicsKind::kSurfaceEnergyBalance) {
              add_surface_budget(pending_surface_budget, step.surface_budget);
              if (derived != nullptr) {
                surface_diagnostics.push_back(
                    {sampled.time_s, step.surface_rates, pending_surface_budget});
                pending_surface_budget = {};
              }
            }
            if (config.physics.kind == mps::PhysicsKind::kGrayRadiation) {
              add_radiation_budget(pending_radiation_budget, step.radiation_budget);
              pending_cfl_retries += step.cfl_retry_count;
              pending_invariant_retries += step.invariant_retry_count;
              pending_solver_retries += step.solver_retry_count;
              pending_radiation_column_calls += step.radiation_column_call_count;
              pending_radiation_wall_seconds += step.radiation_wall_seconds;
              if (derived != nullptr) {
                radiation_diagnostics.push_back(
                    {.time_s = sampled.time_s,
                     .step = sampled.step,
                     .segment_start_time_s = segment_start_time_s,
                     .segment_start_step = segment_start_step,
                     .rates = step.radiation_rates,
                     .interval_budget = pending_radiation_budget,
                     .requested_time_step_s = step.requested_time_step_s,
                     .accepted_time_step_s = step.accepted_time_step_s,
                     .radiation_stable_time_step_s = step.radiation_stable_time_step_s,
                     .cfl_retry_count = pending_cfl_retries,
                     .invariant_retry_count = pending_invariant_retries,
                     .solver_retry_count = pending_solver_retries,
                     .radiation_column_call_count = pending_radiation_column_calls,
                     .radiation_wall_seconds = pending_radiation_wall_seconds});
                pending_radiation_budget = {};
                pending_cfl_retries = 0;
                pending_invariant_retries = 0;
                pending_solver_retries = 0;
                pending_radiation_column_calls = 0;
                pending_radiation_wall_seconds = 0.0;
              }
            }
          },
          [&state, &command_line] {
            return g_cancel_requested.load() ||
                   (command_line.stop_after_step.has_value() &&
                    state.step >= *command_line.stop_after_step);
          },
          [&period_mean, &visual_writer, &visual_fields, &config, &driver](
              const mps::DryHydrostaticState& accepted, const mps::Real step_start_s,
              const mps::Real step_end_s) {
            if (!period_mean.has_value()) return;
            const auto derived = driver.diagnose(accepted);
            const auto values = mps::dry_visual_values(config, driver.grid(), accepted,
                                                       derived, visual_fields);
            period_mean->observe(step_start_s, step_end_s, values);
            for (const auto& completed : period_mean->take_completed())
              visual_writer.publish(completed);
          });
      const bool reached_end_time = state.time_s >= config.run.end_time_s;
      const bool was_cancelled = g_cancel_requested.load();
      const std::string_view result_status =
          reached_end_time ? "complete" : (was_cancelled ? "cancelled" : "stopped");
      progress.finish(state, result_status);
      if (period_mean.has_value()) {
        if (const auto partial = period_mean->partial(); partial.has_value())
          visual_writer.publish(*partial);
      }
      if (physics_diagnostics.has_value()) physics_diagnostics->write();
      if (moist_diagnostics.has_value()) moist_diagnostics->write();
      if (climate_statistics.has_value()) {
        const std::filesystem::path directory(config.output_directory);
        std::filesystem::create_directories(directory);
        std::ofstream output(directory / "climate_statistics.csv", std::ios::trunc);
        if (!output) throw std::runtime_error("unable to open climate statistics CSV");
        mps::write_climate_statistics_csv(output, climate_statistics->rows());
      }
      if (driver.semi_implicit_vertical_modes().has_value())
        write_semi_implicit_diagnostics(config, semi_implicit_diagnostics);
      if (config.convection.kind != mps::ConvectionKind::kNone ||
          config.boundary_layer.kind != mps::BoundaryLayerKind::kNone)
        write_dry_mixing_diagnostics(config, driver.grid(), dry_mixing_diagnostics);
      if (config.physics.kind == mps::PhysicsKind::kSurfaceEnergyBalance ||
          config.physics.kind == mps::PhysicsKind::kGrayRadiation) {
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
        if (config.physics.kind == mps::PhysicsKind::kSurfaceEnergyBalance) {
          std::ofstream diagnostics_output(directory / "surface_diagnostics.csv",
                                           std::ios::trunc);
          if (!diagnostics_output)
            throw std::runtime_error("unable to open surface diagnostics CSV");
          diagnostics_output
              << std::setprecision(17)
              << "time_s,absorbed_stellar_power_w,internal_heat_power_w,"
                 "outgoing_longwave_power_w,sensible_to_atmosphere_power_w,"
                 "surface_storage_rate_w,interval_absorbed_stellar_energy_j,"
                 "interval_internal_heat_energy_j,"
                 "interval_outgoing_longwave_energy_j,"
                 "interval_sensible_to_atmosphere_energy_j,"
                 "interval_surface_storage_change_j,"
                 "interval_surface_budget_residual_j\n";
          for (const auto& sample : surface_diagnostics) {
            const auto& rates = sample.rates;
            const auto& budget = sample.interval_budget;
            diagnostics_output << sample.time_s << ',' << rates.absorbed_stellar_power_w
                               << ',' << rates.internal_heat_power_w << ','
                               << rates.outgoing_longwave_power_w << ','
                               << rates.sensible_to_atmosphere_power_w << ','
                               << rates.surface_storage_rate_w << ','
                               << budget.absorbed_stellar_energy_j << ','
                               << budget.internal_heat_energy_j << ','
                               << budget.outgoing_longwave_energy_j << ','
                               << budget.sensible_to_atmosphere_energy_j << ','
                               << budget.surface_storage_change_j << ','
                               << budget.surface_budget_residual_j << '\n';
          }
        }
      }
      if (config.physics.kind == mps::PhysicsKind::kGrayRadiation)
        write_radiation_diagnostics(config, driver.grid(), radiation_diagnostics);
      if (command_line.checkpoint_path.has_value()) {
        const bool has_surface =
            config.physics.kind == mps::PhysicsKind::kSurfaceEnergyBalance ||
            config.physics.kind == mps::PhysicsKind::kGrayRadiation;
        const bool multitracer = !config.tracers.empty();
        const bool moist = config.moisture.kind == mps::MoistureKind::kDiluteWater;
        const mps::TracerRegistry registry = multitracer
                                                 ? mps::TracerRegistry(config.tracers)
                                                 : mps::TracerRegistry::legacy();
        mps::write_checkpoint_file(
            *command_line.checkpoint_path,
            {.time_s = state.time_s,
             .step = state.step,
             .state =
                 moist ? mps::flatten_dry_hydrostatic_moist_state(
                             state, static_cast<std::size_t>(config.vertical.levels))
                 : multitracer
                     ? mps::flatten_dry_hydrostatic_multitracer_state(
                           state, static_cast<std::size_t>(config.vertical.levels),
                           has_surface)
                 : has_surface
                     ? mps::flatten_dry_hydrostatic_surface_state(
                           state, static_cast<std::size_t>(config.vertical.levels))
                     : mps::flatten_dry_hydrostatic_state(
                           state, static_cast<std::size_t>(config.vertical.levels)),
             .config_fingerprint = fingerprint,
             .layout_id =
                 moist ? mps::dry_hydrostatic_moist_checkpoint_layout(registry)
                 : multitracer
                     ? mps::dry_hydrostatic_multitracer_checkpoint_layout(registry,
                                                                          has_surface)
                     : std::string(has_surface
                                       ? mps::kDryHydrostaticSurfaceCheckpointLayout
                                       : mps::kDryHydrostaticCheckpointLayout)});
        if (period_mean.has_value()) {
          mps::write_statistics_sidecar_file(
              mps::statistics_sidecar_path(*command_line.checkpoint_path),
              {.checkpoint_hash = mps::file_fnv1a64(*command_line.checkpoint_path),
               .checkpoint_time_s = state.time_s,
               .checkpoint_step = state.step,
               .accumulator = period_mean->state(),
               .periods = visual_writer.periods()});
        }
      }
      const auto derived = driver.diagnose(state);
      if (config.physics.kind == mps::PhysicsKind::kGrayRadiation)
        write_radiation_column(config, driver, state, derived);
      if (config.boundary_layer.kind != mps::BoundaryLayerKind::kNone)
        write_mixing_column(config, driver, state, derived);
      if (config.moisture.kind == mps::MoistureKind::kDiluteWater)
        write_moist_column(config, driver, state, derived);
      const auto diagnostics = mps::diagnose_dry_hydrostatic_budgets(
          driver.grid(), state, derived, config.planet);
      mps::write_run_metadata(std::cout, mps::make_run_metadata(config), config);
      if (driver.semi_implicit_vertical_modes().has_value())
        mps::write_dry_hydrostatic_vertical_mode_metadata(
            std::cout, *driver.semi_implicit_vertical_modes());
      if (config.orography.kind != mps::OrographyKind::kFlat) {
        const auto sources = driver.diagnose_sources(derived);
        const auto terrain = mps::diagnose_terrain_budgets(
            driver.grid(), state, derived, sources,
            driver.orography().surface_geopotential_m2_s2(), config.planet);
        mps::write_terrain_diagnostics(std::cout, terrain);
      }
      std::cout << "result.status = " << result_status
                << "\nresult.time_s = " << state.time_s
                << "\nresult.step = " << state.step
                << "\ndiagnostics.dry_mass_kg = " << diagnostics.dry_mass_kg
                << "\ndiagnostics.total_energy_j = " << diagnostics.total_energy_j
                << '\n';
      if (was_cancelled) return 130;
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
      const auto start = std::chrono::steady_clock::now();
      const auto result = mps::run_shallow_water(
          config, std::move(initial_state), command_line.stop_after_step);
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
      write_shallow_water_snapshot(config, result);
      write_shallow_water_diagnostics(config, result);
      mps::write_run_metadata(std::cout, mps::make_run_metadata(config), config);
      write_shallow_water_result(std::cout, result, wall_time_s, grid.cell_count());
      write_shallow_water_errors(std::cout, config, grid, result);
    }
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << '\n';
    return 1;
  }
}
