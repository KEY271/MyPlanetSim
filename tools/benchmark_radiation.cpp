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
    config.diagnostics.interval_steps = static_cast<std::uint64_t>(
        std::max(1.0, std::ceil(86400.0 / config.run.time_step_s)));
    if (iterations > 0) {
      const auto solver_parameters =
          config.semi_implicit.value_or(mps::SemiImplicitParameters{});
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
          .linear_relative_tolerance =
              config.semi_implicit ? solver_parameters.linear_relative_tolerance : 1e-8,
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
    std::ofstream step_output(std::string(argv[8]) + ".steps.csv");
    std::ofstream climate_output(std::string(argv[8]) + ".climate.csv");
    if (!step_output || !climate_output)
      throw std::runtime_error("cannot open pilot diagnostics");
    step_output
        << "day,step,dt_s,linear_iterations,cfl_retries,invariant_retries,solver_"
           "retries,rhs_s,linear_s,total_s,toa_net_j,surface_storage_j,dry_thermal_j,"
           "drag_j,diffusion_j,interface_residual_j,surface_integration_residual_j\n"
        << std::setprecision(17);
    climate_output << "day,mean_temperature_k,mean_surface_temperature_k,maximum_wind_"
                      "m_s,toa_net_w_m2,surface_storage_w_m2,unstable_fraction,minimum_"
                      "dtheta_dz_k_m,mass_relative_drift,energy_residual_j\n"
                   << std::setprecision(17);
    double attributed_energy = 0.0;
    double maximum_mass_drift = 0.0;
    const auto start = std::chrono::steady_clock::now();
    try {
      driver.advance(
          state, config.run.end_time_s,
          [&](const mps::DryHydrostaticState& sampled,
              const mps::DryHydrostaticDerived* sampled_derived,
              const mps::DryHydrostaticStepDiagnostics& step) {
            if (sampled_derived != nullptr) {
              const auto& d = *sampled_derived;
              double mean_temperature = 0.0, mean_surface = 0.0, wind = 0.0;
              double unstable = 0.0, minimum_gradient = 0.0;
              for (std::size_t cell = 0; cell < driver.grid().cell_count(); ++cell) {
                const auto area = driver.grid().cells()[cell].area_m2;
                mean_surface += area * sampled.surface_temperature_k[cell];
                for (int k = 0; k < levels; ++k) {
                  const auto n = cell * levels + k;
                  mean_temperature += area * d.temperature_k[n] / levels;
                  wind = std::max(wind, mps::norm(d.velocity_m_s[n]));
                  if (k + 1 < levels) {
                    const auto dz =
                        (d.geopotential_m2_s2[n] - d.geopotential_m2_s2[n + 1]) /
                        config.planet.gravity_m_s2;
                    const auto gradient = (d.potential_temperature_k[n] -
                                           d.potential_temperature_k[n + 1]) /
                                          dz;
                    minimum_gradient = std::min(minimum_gradient, gradient);
                    if (gradient < 0.0) unstable += area / (levels - 1);
                  }
                }
              }
              const auto health = mps::diagnose_dry_hydrostatic_budgets(
                  driver.grid(), sampled, d, config.planet);
              const double drift =
                  (health.dry_mass_kg - initial.dry_mass_kg) / initial.dry_mass_kg;
              maximum_mass_drift = std::max(maximum_mass_drift, std::abs(drift));
              const auto area = driver.grid().total_area_m2();
              const auto current_energy = step.thermal_energy_contribution_j +
                                          step.rayleigh_drag_energy_contribution_j +
                                          step.diffusion_energy_contribution_j;
              climate_output << sampled.time_s / 86400.0 << ','
                             << mean_temperature / area << ',' << mean_surface / area
                             << ',' << wind << ','
                             << step.radiation_rates.toa_net_upward_power_w / area
                             << ','
                             << step.radiation_rates.surface_storage_rate_w / area
                             << ',' << unstable / area << ',' << minimum_gradient << ','
                             << drift << ','
                             << health.total_energy_j - initial.total_energy_j -
                                    attributed_energy - current_energy
                             << '\n';
            }
            if (step.accepted_time_step_s == 0.0) return;
            attributed_energy += step.thermal_energy_contribution_j +
                                 step.rayleigh_drag_energy_contribution_j +
                                 step.diffusion_energy_contribution_j;
            const auto& budget = step.radiation_budget;
            step_output << sampled.time_s / 86400.0 << ',' << sampled.step << ','
                        << step.accepted_time_step_s << ','
                        << step.linear_iterations_total << ',' << step.cfl_retry_count
                        << ',' << step.invariant_retry_count << ','
                        << step.solver_retry_count << ',' << step.wall_seconds_rhs
                        << ',' << step.wall_seconds_linear_solve << ','
                        << step.wall_seconds_total << ','
                        << budget.toa_net_upward_energy_j << ','
                        << budget.surface_storage_change_j << ','
                        << budget.dry_thermal_energy_j << ','
                        << budget.rayleigh_drag_energy_j << ','
                        << step.diffusion_energy_contribution_j << ','
                        << budget.interface_conservation_residual_j << ','
                        << budget.surface_time_integration_residual_j << '\n';
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
    if (!step_output || !climate_output)
      throw std::runtime_error("cannot write pilot diagnostics");
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
              << "maximum_sampled_mass_drift=" << maximum_mass_drift << '\n'
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
