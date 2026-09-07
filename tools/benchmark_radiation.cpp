#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

#if defined(__unix__) || defined(__APPLE__)
#include <sys/resource.h>
#endif

#include "myplanetsim/diagnostics/dry_hydrostatic_diagnostics.hpp"
#include "myplanetsim/dynamics/dry_hydrostatic_driver.hpp"

// One process per experiment keeps timing independent of other matrix members.
// Snapshot rows retain every atmospheric and surface temperature for RMS comparisons.
int main(int argc, char** argv) {
  if (argc != 9) {
    std::cerr
        << "usage: benchmark_radiation CONFIG DAYS N K DT ITERATIONS TOP_PA SNAPSHOT\n"
           "ITERATIONS=0 selects SSP-RK3; positive values select centered "
           "semi-implicit.\n";
    return 2;
  }
  try {
    auto config = mps::load_experiment_config(argv[1]);
    config.run.end_time_s = config.run.start_time_s + std::stod(argv[2]) * 86400.0;
    config.grid.cells_per_panel = std::stoi(argv[3]);
    const int levels = std::stoi(argv[4]);
    if (levels < 1 || levels > 1000) throw std::invalid_argument("invalid level count");
    config.run.time_step_s = std::stod(argv[5]);
    const int iterations = std::stoi(argv[6]);
    if (iterations < 0) throw std::invalid_argument("negative iteration count");
    const double top = std::stod(argv[7]);
    config.vertical.levels = levels;
    config.vertical.a_half_pa.resize(levels + 1);
    config.vertical.b_half.resize(levels + 1);
    for (int k = 0; k <= levels; ++k) {
      const double fraction = static_cast<double>(k) / levels;
      config.vertical.a_half_pa[k] = top * (1.0 - fraction);
      config.vertical.b_half[k] = fraction;
    }
    config.diagnostics.interval_steps = 1000000;
    if (iterations > 0) {
      config.dry_hydrostatic.time_integrator =
          mps::DryHydrostaticTimeIntegrator::kSemiImplicit;
      config.dry_hydrostatic.advective_cfl = 0.45;
      config.semi_implicit = mps::SemiImplicitParameters{
          .reference_surface_pressure_pa = config.vertical.surface_pressure_pa,
          .reference_temperature_k = config.vertical.initial_temperature_k,
          .reference_update = mps::SemiImplicitReferenceUpdate::kFixed,
          .implicit_weight = 0.5,
          .wave_cfl_threshold = 0.45,
          .maximum_implicit_modes = std::min(levels, 5),
          .nonlinear_iterations = iterations,
          .linear_relative_tolerance = 1e-8,
          .linear_absolute_tolerance = 1e-12,
          .linear_maximum_iterations = 80,
          .gmres_restart = 20,
          .minimum_time_step_s = 0.01};
    } else {
      config.dry_hydrostatic.time_integrator =
          mps::DryHydrostaticTimeIntegrator::kExplicitSspRk3;
      config.semi_implicit.reset();
    }
    config.validate();
    const mps::DryHydrostaticDriver driver(config);
    auto state = driver.initial_state();
    const auto initial = mps::diagnose_dry_hydrostatic_budgets(
        driver.grid(), state, driver.diagnose(state), config.planet);
    std::vector<double> steps;
    double radiation_seconds = 0.0;
    std::size_t calls = 0, retries = 0, linear_iterations = 0;
    double max_interface_residual = 0.0;
    const auto start = std::chrono::steady_clock::now();
    try {
      driver.advance(
          state, config.run.end_time_s,
          [&](const mps::DryHydrostaticState& sampled,
              const mps::DryHydrostaticDerived*,
              const mps::DryHydrostaticStepDiagnostics& step) {
            if (step.accepted_time_step_s == 0.0) return;
            steps.push_back(step.accepted_time_step_s);
            radiation_seconds += step.radiation_wall_seconds;
            calls += step.radiation_column_call_count;
            retries += step.retry_count;
            linear_iterations += step.linear_iterations_total;
            max_interface_residual = std::max(
                max_interface_residual,
                std::abs(step.radiation_budget.interface_conservation_residual_j) /
                    (step.accepted_time_step_s * driver.grid().total_area_m2()));
            if (sampled.step % 1000 == 0)
              std::cerr << "day=" << sampled.time_s / 86400.0
                        << " step=" << sampled.step << '\n';
          });
    } catch (const std::exception& error) {
      std::cerr << "failed_day=" << state.time_s / 86400.0 << " step=" << state.step
                << " reason=" << error.what() << '\n';
      throw;
    }
    const double elapsed =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    if (steps.empty()) throw std::runtime_error("no accepted steps");
    const auto derived = driver.diagnose(state);
    const auto final = mps::diagnose_dry_hydrostatic_budgets(driver.grid(), state,
                                                             derived, config.planet);
    mps::DryHydrostaticRhs rhs;
    driver.rhs(state, rhs);
    std::sort(steps.begin(), steps.end());
    std::ofstream snapshot(argv[8]);
    if (!snapshot) throw std::runtime_error("cannot open snapshot");
    snapshot << "cell,level,area_m2,temperature_k,surface_temperature_k\n"
             << std::setprecision(17);
    for (std::size_t cell = 0; cell < driver.grid().cell_count(); ++cell)
      for (int k = 0; k < levels; ++k)
        snapshot << cell << ',' << k << ',' << driver.grid().cells()[cell].area_m2
                 << ',' << derived.temperature_k[cell * levels + k] << ','
                 << state.surface_temperature_k[cell] << '\n';
    if (!snapshot) throw std::runtime_error("cannot write snapshot");
    std::size_t peak_rss_bytes = 0;
#if defined(__unix__) || defined(__APPLE__)
    rusage usage{};
    if (getrusage(RUSAGE_SELF, &usage) == 0) {
      peak_rss_bytes = static_cast<std::size_t>(usage.ru_maxrss);
#if !defined(__APPLE__)
      peak_rss_bytes *= 1024;
#endif
    }
#endif
    std::cout << std::setprecision(17) << "peak_rss_bytes=" << peak_rss_bytes << '\n'
              << "elapsed_s=" << elapsed << '\n'
              << "seconds_per_model_day="
              << elapsed * 86400.0 / (state.time_s - config.run.start_time_s) << '\n'
              << "accepted_steps=" << steps.size() << '\n'
              << "median_dt_s=" << steps[steps.size() / 2] << '\n'
              << "minimum_dt_s=" << steps.front() << '\n'
              << "radiation_seconds=" << radiation_seconds << '\n'
              << "radiation_calls=" << calls << '\n'
              << "retries=" << retries << '\n'
              << "linear_iterations=" << linear_iterations << '\n'
              << "mass_relative_drift="
              << (final.dry_mass_kg - initial.dry_mass_kg) / initial.dry_mass_kg << '\n'
              << "interface_residual_w_m2=" << max_interface_residual << '\n'
              << "toa_net_upward_w_m2="
              << rhs.radiation_diagnostics.toa_net_upward_power_w /
                     driver.grid().total_area_m2()
              << '\n'
              << "minimum_temperature_k=" << final.minimum_temperature_k << '\n'
              << "maximum_temperature_k=" << final.maximum_temperature_k << '\n';
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
