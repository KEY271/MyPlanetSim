#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#if defined(__unix__) || defined(__APPLE__)
#include <sys/resource.h>
#endif

#include "myplanetsim/diagnostics/dry_hydrostatic_diagnostics.hpp"
#include "myplanetsim/dynamics/dry_hydrostatic_driver.hpp"
#include "myplanetsim/vertical/hybrid_pressure_coordinate.hpp"

// One process per experiment, as in benchmark_radiation, so timings stay independent.
// The vertical coordinate is the surface-refined ramp Phase 12 uses; `REFINEMENT=1`
// reproduces the uniform sigma grid of the earlier phases for a controlled comparison.
// LEG selects the physics combination of the plan's four-way matrix. Legs without a
// boundary layer keep the Phase 11 owner of the sensible heat flux and Rayleigh drag.
int main(int argc, char** argv) {
  if (argc != 11) {
    std::cerr << "usage: benchmark_dry_mixing CONFIG DAYS N K DT ITERATIONS TOP_PA "
                 "REFINEMENT LEG SNAPSHOT\n"
                 "ITERATIONS=0 selects SSP-RK3; positive values select centered "
                 "semi-implicit.\n"
                 "LEG is one of radiation, convection, boundary_layer, all.\n";
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
    const double refinement = std::stod(argv[8]);
    const std::string leg = argv[9];
    const bool boundary_layer = leg == "boundary_layer" || leg == "all";
    const bool convection = leg == "convection" || leg == "all";
    if (!boundary_layer && !convection && leg != "radiation")
      throw std::invalid_argument("unknown physics leg " + leg);

    const auto coefficients = mps::surface_refined_sigma_coefficients(top, levels,
                                                                     refinement);
    config.vertical.levels = levels;
    config.vertical.a_half_pa = coefficients.a_half_pa;
    config.vertical.b_half = coefficients.b_half;
    config.boundary_layer.kind = boundary_layer ? mps::BoundaryLayerKind::kBulkKProfile
                                                : mps::BoundaryLayerKind::kNone;
    config.convection.kind =
        convection ? mps::ConvectionKind::kDryAdjustment : mps::ConvectionKind::kNone;
    config.surface->air_exchange_coefficient_w_m2_k = boundary_layer ? 0.0 : 10.0;
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
              config.semi_implicit ? solver_parameters.linear_relative_tolerance : 1e-10,
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
    const double area = driver.grid().total_area_m2();
    auto state = driver.initial_state();
    const auto initial = mps::diagnose_dry_hydrostatic_budgets(
        driver.grid(), state, driver.diagnose(state), config.planet);
    std::vector<double> steps;
    std::size_t retries = 0, linear_iterations = 0;
    std::size_t boundary_layer_calls = 0, convection_calls = 0;
    double maximum_heat_residual = 0.0, maximum_kinetic_residual = 0.0;
    double maximum_tracer_change = 0.0, maximum_momentum_residual = 0.0;
    double maximum_unstable_before = 0.0, maximum_unstable_after = 0.0;
    double maximum_boundary_layer_height = 0.0;
    double sensible_energy = 0.0, dissipation_energy = 0.0;
    double surface_boundary_layer_energy = 0.0;
    std::ofstream step_output(std::string(argv[10]) + ".steps.csv");
    std::ofstream climate_output(std::string(argv[10]) + ".climate.csv");
    if (!step_output || !climate_output)
      throw std::runtime_error("cannot open pilot diagnostics");
    step_output << "day,step,dt_s,linear_iterations,cfl_retries,invariant_retries,"
                   "solver_retries,rhs_s,total_s,toa_net_j,surface_storage_j,"
                   "legacy_sensible_j,legacy_drag_j,boundary_layer_sensible_j,"
                   "boundary_layer_surface_j,returned_dissipation_j,"
                   "heat_residual_j,kinetic_identity_residual_j,"
                   "momentum_residual_n_s,tracer_change_kg,convection_enthalpy_j,"
                   "convection_dry_energy_j,unstable_before,unstable_after,mean_h_m,"
                   "maximum_k_h_m2_s,shallow_fraction,model_top_fraction\n"
                << std::setprecision(17);
    climate_output << "day,mean_temperature_k,mean_surface_temperature_k,maximum_wind_"
                      "m_s,toa_net_w_m2,surface_storage_w_m2,unstable_fraction,minimum_"
                      "dtheta_dz_k_m,mean_h_m,mass_relative_drift,unattributed_energy_j\n"
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
                const auto cell_area = driver.grid().cells()[cell].area_m2;
                mean_surface += cell_area * sampled.surface_temperature_k[cell];
                for (int k = 0; k < levels; ++k) {
                  const auto n = cell * static_cast<std::size_t>(levels) + k;
                  mean_temperature += cell_area * d.temperature_k[n] / levels;
                  wind = std::max(wind, mps::norm(d.velocity_m_s[n]));
                  if (k + 1 < levels) {
                    const auto dz =
                        (d.geopotential_m2_s2[n] - d.geopotential_m2_s2[n + 1]) /
                        config.planet.gravity_m_s2;
                    const auto gradient = (d.potential_temperature_k[n] -
                                           d.potential_temperature_k[n + 1]) /
                                          dz;
                    minimum_gradient = std::min(minimum_gradient, gradient);
                    if (gradient < 0.0) unstable += cell_area / (levels - 1);
                  }
                }
              }
              const auto health = mps::diagnose_dry_hydrostatic_budgets(
                  driver.grid(), sampled, d, config.planet);
              const double drift =
                  (health.dry_mass_kg - initial.dry_mass_kg) / initial.dry_mass_kg;
              maximum_mass_drift = std::max(maximum_mass_drift, std::abs(drift));
              const auto current_energy = step.thermal_energy_contribution_j +
                                          step.rayleigh_drag_energy_contribution_j +
                                          step.diffusion_energy_contribution_j;
              climate_output << sampled.time_s / 86400.0 << ','
                             << mean_temperature / area << ',' << mean_surface / area
                             << ',' << wind << ','
                             << step.radiation_rates.toa_net_upward_power_w / area << ','
                             << step.radiation_rates.surface_storage_rate_w / area << ','
                             << unstable / area << ',' << minimum_gradient << ','
                             << step.boundary_layer.mean_boundary_layer_height_m << ','
                             << drift << ','
                             << health.total_energy_j - initial.total_energy_j -
                                    attributed_energy - current_energy
                             << '\n';
            }
            if (step.accepted_time_step_s == 0.0) return;
            // The convective adjustment reports an exact dry-energy attribution, so it
            // is subtracted here. The boundary layer reports its change on the enthalpy
            // basis rather than the cv+Phi basis of `total_energy_j`, so its effect is
            // deliberately left inside the residual column; that column is therefore
            // `unattributed`, not a defect measurement, and is not comparable to the
            // Phase 11 residual.
            attributed_energy += step.thermal_energy_contribution_j +
                                 step.rayleigh_drag_energy_contribution_j +
                                 step.diffusion_energy_contribution_j +
                                 step.convection.dry_energy_attributed_change_j;
            const auto& budget = step.radiation_budget;
            const auto& mixing = step.boundary_layer;
            step_output << sampled.time_s / 86400.0 << ',' << sampled.step << ','
                        << step.accepted_time_step_s << ','
                        << step.linear_iterations_total << ',' << step.cfl_retry_count
                        << ',' << step.invariant_retry_count << ','
                        << step.solver_retry_count << ',' << step.wall_seconds_rhs << ','
                        << step.wall_seconds_total << ','
                        << budget.toa_net_upward_energy_j << ','
                        << budget.surface_storage_change_j << ','
                        << budget.sensible_to_atmosphere_energy_j << ','
                        << budget.rayleigh_drag_energy_j << ','
                        << mixing.sensible_to_atmosphere_energy_j << ','
                        << mixing.surface_heat_change_j << ','
                        << mixing.returned_dissipation_heat_j << ','
                        << mixing.heat_budget_residual_j << ','
                        << mixing.kinetic_energy_identity_residual_j << ','
                        << mixing.momentum_budget_residual_n_s << ','
                        << mixing.tracer_mass_change_kg << ','
                        << step.convection.enthalpy_change_j << ','
                        << step.convection.dry_energy_attributed_change_j << ','
                        << step.convection.unstable_interface_fraction_before << ','
                        << step.convection.unstable_interface_fraction_after << ','
                        << mixing.mean_boundary_layer_height_m << ','
                        << mixing.maximum_heat_diffusivity_m2_s << ','
                        << mixing.shallow_unresolved_area_fraction << ','
                        << mixing.model_top_area_fraction << '\n';
            steps.push_back(step.accepted_time_step_s);
            retries += step.retry_count;
            linear_iterations += step.linear_iterations_total;
            boundary_layer_calls += step.boundary_layer_column_call_count;
            convection_calls += step.convection_column_call_count;
            sensible_energy += mixing.sensible_to_atmosphere_energy_j;
            dissipation_energy += mixing.returned_dissipation_heat_j;
            surface_boundary_layer_energy += mixing.surface_heat_change_j;
            const double heat_scale =
                std::max(1.0 * area, std::abs(mixing.atmospheric_heat_change_j) +
                                         std::abs(mixing.surface_heat_change_j) +
                                         std::abs(mixing.returned_dissipation_heat_j));
            maximum_heat_residual = std::max(
                maximum_heat_residual, std::abs(mixing.heat_budget_residual_j) /
                                           heat_scale);
            maximum_kinetic_residual = std::max(
                maximum_kinetic_residual,
                std::abs(mixing.kinetic_energy_identity_residual_j) /
                    std::max(1.0 * area, std::abs(mixing.kinetic_energy_change_j)));
            maximum_tracer_change =
                std::max(maximum_tracer_change, std::abs(mixing.tracer_mass_change_kg));
            maximum_momentum_residual =
                std::max(maximum_momentum_residual, mixing.momentum_budget_residual_n_s);
            maximum_unstable_before =
                std::max(maximum_unstable_before,
                         step.convection.unstable_interface_fraction_before);
            maximum_unstable_after =
                std::max(maximum_unstable_after,
                         step.convection.unstable_interface_fraction_after);
            maximum_boundary_layer_height = std::max(
                maximum_boundary_layer_height, mixing.mean_boundary_layer_height_m);
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
    std::sort(steps.begin(), steps.end());
    std::ofstream snapshot(argv[10]);
    if (!snapshot) throw std::runtime_error("cannot open snapshot");
    snapshot << "cell,level,area_m2,air_mass_kg_m2,temperature_k,"
                "potential_temperature_k,wind_m_s,surface_temperature_k\n"
             << std::setprecision(17);
    for (std::size_t cell = 0; cell < driver.grid().cell_count(); ++cell)
      for (int k = 0; k < levels; ++k) {
        const auto n = cell * static_cast<std::size_t>(levels) + k;
        snapshot << cell << ',' << k << ',' << driver.grid().cells()[cell].area_m2 << ','
                 << derived.air_mass_kg_m2[n] << ',' << derived.temperature_k[n] << ','
                 << derived.potential_temperature_k[n] << ','
                 << mps::norm(derived.velocity_m_s[n]) << ','
                 << state.surface_temperature_k[cell] << '\n';
      }
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
              << "retries=" << retries << '\n'
              << "linear_iterations=" << linear_iterations << '\n'
              << "boundary_layer_column_calls=" << boundary_layer_calls << '\n'
              << "convection_column_calls=" << convection_calls << '\n'
              << "sensible_energy_j=" << sensible_energy << '\n'
              << "boundary_layer_surface_energy_j=" << surface_boundary_layer_energy
              << '\n'
              << "returned_dissipation_j=" << dissipation_energy << '\n'
              << "maximum_heat_residual=" << maximum_heat_residual << '\n'
              << "maximum_kinetic_identity_residual=" << maximum_kinetic_residual << '\n'
              << "maximum_tracer_change_kg=" << maximum_tracer_change << '\n'
              << "maximum_momentum_residual_n_s=" << maximum_momentum_residual << '\n'
              << "maximum_unstable_before=" << maximum_unstable_before << '\n'
              << "maximum_unstable_after=" << maximum_unstable_after << '\n'
              << "maximum_mean_boundary_layer_height_m="
              << maximum_boundary_layer_height << '\n'
              << "maximum_sampled_mass_drift=" << maximum_mass_drift << '\n'
              << "mass_relative_drift="
              << (final.dry_mass_kg - initial.dry_mass_kg) / initial.dry_mass_kg << '\n'
              << "minimum_temperature_k=" << final.minimum_temperature_k << '\n'
              << "maximum_temperature_k=" << final.maximum_temperature_k << '\n';
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
