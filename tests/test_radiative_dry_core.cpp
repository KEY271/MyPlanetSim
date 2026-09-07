#include <algorithm>
#include <cmath>
#include <vector>

#include "myplanetsim/diagnostics/dry_hydrostatic_diagnostics.hpp"
#include "myplanetsim/dynamics/dry_hydrostatic_driver.hpp"
#include "support/test.hpp"

namespace {

[[nodiscard]] mps::ExperimentConfig gray_config() {
  mps::ExperimentConfig config{
      .kind = mps::ExperimentKind::kDryHydrostatic,
      .planet = {6.37122e6, 7.292115e-5, 9.80616, 287.0, 1004.0, 100000.0},
      .run = {0.0, 60.0, 10.0, 7},
      .grid = {2},
      .vertical = {.levels = 2,
                   .a_half_pa = {1000.0, 500.0, 0.0},
                   .b_half = {0.0, 0.5, 1.0},
                   .surface_pressure_pa = 100000.0,
                   .minimum_surface_pressure_pa = 90000.0,
                   .maximum_surface_pressure_pa = 110000.0,
                   .minimum_pressure_thickness_pa = 100.0,
                   .initial_temperature_k = 280.0,
                   .initial_potential_temperature_k = 300.0,
                   .temperature_floor_k = 100.0,
                   .transport_scheme = mps::VerticalTransportScheme::kDonorCell,
                   .limiter = mps::VerticalLimiterKind::kNone,
                   .cfl = 0.5},
      .dry_hydrostatic = {},
      .physics = {.kind = mps::PhysicsKind::kGrayRadiation},
      .diagnostics = {1},
      .output_directory = "x",
      .orbit =
          mps::OrbitParameters{
              .period_s = 86400.0,
              .eccentricity = 0.0,
              .obliquity_rad = 0.0,
              .longitude_of_periapsis_rad = 0.0,
              .initial_mean_anomaly_rad = 0.0,
              .initial_substellar_longitude_rad = 0.0,
              .stellar_flux_at_semimajor_axis_w_m2 = 1000.0,
          },
      .surface =
          mps::SurfaceParameters{
              .geography = mps::SurfaceGeography::kUniform,
              .uniform_land_fraction = 0.25,
              .land_heat_capacity_j_m2_k = 2e6,
              .ocean_heat_capacity_j_m2_k = 2e8,
              .initial_temperature_k = 290.0,
              .albedo = 0.25,
              .emissivity = 0.95,
              .air_exchange_coefficient_w_m2_k = 8.0,
              .internal_heat_flux_w_m2 = 0.1,
              .cfl = 0.5,
          },
      .radiation =
          mps::RadiationParameters{
              .shortwave_absorption_m2_kg = 2e-5,
              .longwave_absorption_ref_m2_kg = 3e-4,
              .reference_pressure_pa = 100000.0,
              .longwave_pressure_exponent = 1.0,
              .longwave_diffusivity_factor = 1.66,
              .shortwave_diffuse_factor = 1.66,
              .cfl = 0.5,
          },
  };
  config.validate();
  return config;
}

[[nodiscard]] mps::Real dry_mass(const mps::DryHydrostaticDriver& driver,
                                 const mps::DryHydrostaticState& state) {
  const auto derived = driver.diagnose(state);
  mps::Real result = 0.0;
  for (std::size_t cell = 0; cell < driver.grid().cell_count(); ++cell)
    for (std::size_t level = 0; level < derived.levels; ++level)
      result +=
          driver.grid().cells()[cell].area_m2 *
          derived
              .air_mass_kg_m2[mps::dry_hydrostatic_offset(cell, level, derived.levels)];
  return result;
}

[[nodiscard]] mps::ExperimentConfig semi_implicit_config(const mps::Real time_step_s,
                                                         const mps::Index iterations) {
  auto config = gray_config();
  config.run.end_time_s = 80.0;
  config.run.time_step_s = time_step_s;
  config.dry_hydrostatic.time_integrator =
      mps::DryHydrostaticTimeIntegrator::kSemiImplicit;
  config.dry_hydrostatic.advective_cfl = 0.45;
  config.semi_implicit = mps::SemiImplicitParameters{
      .reference_surface_pressure_pa = config.vertical.surface_pressure_pa,
      .reference_temperature_k = config.vertical.initial_temperature_k,
      .reference_update = mps::SemiImplicitReferenceUpdate::kFixed,
      .implicit_weight = 0.5,
      .wave_cfl_threshold = 0.45,
      .maximum_implicit_modes = 2,
      .nonlinear_iterations = iterations,
      .linear_relative_tolerance = 1e-8,
      .linear_absolute_tolerance = 1e-12,
      .linear_maximum_iterations = 80,
      .gmres_restart = 20,
      .minimum_time_step_s = 0.01,
  };
  config.validate();
  return config;
}

[[nodiscard]] mps::ExperimentConfig convection_config(const bool semi_implicit) {
  auto config = semi_implicit ? semi_implicit_config(10.0, 5) : gray_config();
  config.run.end_time_s = 20.0;
  config.run.time_step_s = 10.0;
  config.convection.kind = mps::ConvectionKind::kDryAdjustment;
  config.convection.stability_tolerance_k = 1e-10;
  config.validate();
  return config;
}

void impose_unstable_columns(const mps::DryHydrostaticDriver& driver,
                             mps::DryHydrostaticState& state) {
  const auto derived = driver.diagnose(state);
  for (std::size_t cell = 0; cell < derived.cells; ++cell) {
    const auto upper = mps::dry_hydrostatic_offset(cell, 0, derived.levels);
    const auto lower = mps::dry_hydrostatic_offset(cell, 1, derived.levels);
    state.potential_temperature_mass_k_kg_m2[upper] =
        derived.air_mass_kg_m2[upper] * 280.0;
    state.potential_temperature_mass_k_kg_m2[lower] =
        derived.air_mass_kg_m2[lower] * 340.0;
  }
}

[[nodiscard]] mps::DryHydrostaticState run(const mps::ExperimentConfig& config) {
  const mps::DryHydrostaticDriver driver(config);
  auto state = driver.initial_state();
  driver.advance(state, config.run.end_time_s);
  return state;
}

[[nodiscard]] mps::Real normalized_difference(const mps::DryHydrostaticState& first,
                                              const mps::DryHydrostaticState& second,
                                              const std::size_t levels) {
  const auto first_values = mps::flatten_dry_hydrostatic_surface_state(first, levels);
  const auto second_values = mps::flatten_dry_hydrostatic_surface_state(second, levels);
  mps::Real squared = 0.0;
  for (std::size_t index = 0; index < first_values.size(); ++index) {
    const mps::Real scale = std::max(1.0, std::abs(first_values[index]));
    const mps::Real difference = (first_values[index] - second_values[index]) / scale;
    squared += difference * difference;
  }
  return std::sqrt(squared / static_cast<mps::Real>(first_values.size()));
}

}  // namespace

MPS_TEST_CASE("gray radiation RHS changes heat but not mass or tracer") {
  const auto config = gray_config();
  const mps::DryHydrostaticDriver driver(config);
  const auto state = driver.initial_state();
  const auto rhs = driver.rhs(state);
  const auto components = driver.rhs_components(state);
  MPS_CHECK_EQ(rhs.surface_temperature_k_s.size(), driver.grid().cell_count());
  MPS_CHECK(std::isfinite(rhs.radiation_stable_time_step_s));
  MPS_CHECK(rhs.radiation_stable_time_step_s > 0.0);
  MPS_CHECK(std::abs(rhs.radiation_diagnostics.interface_conservation_residual_w) <
            1e-12 * config.orbit->stellar_flux_at_semimajor_axis_w_m2 *
                driver.grid().total_area_m2());
  for (const mps::Real rate : components.physics.surface_pressure_pa_s)
    MPS_CHECK_NEAR(rate, 0.0, 0.0);
  for (const mps::Real rate : components.physics.tendency.air_mass)
    MPS_CHECK_NEAR(rate, 0.0, 0.0);
  for (const mps::Real rate : components.physics.tendency.tracer_mass)
    MPS_CHECK_NEAR(rate, 0.0, 0.0);
  MPS_CHECK(std::ranges::any_of(components.physics.tendency.potential_temperature_mass,
                                [](const mps::Real rate) { return rate != 0.0; }));
}

MPS_TEST_CASE("gray thermal energy attribution matches diagnosed direction") {
  const auto config = gray_config();
  const mps::DryHydrostaticDriver driver(config);
  const auto state = driver.initial_state();
  const auto components = driver.rhs_components(state);
  const auto rhs = driver.rhs(state);
  constexpr mps::Real epsilon = 1e-3;
  auto upper = state;
  auto lower = state;
  for (std::size_t index = 0;
       index < components.physics.tendency.potential_temperature_mass.size(); ++index) {
    const mps::Real change =
        epsilon * components.physics.tendency.potential_temperature_mass[index];
    upper.potential_temperature_mass_k_kg_m2[index] += change;
    lower.potential_temperature_mass_k_kg_m2[index] -= change;
  }
  const auto upper_energy = mps::diagnose_dry_hydrostatic_budgets(
      driver.grid(), upper, driver.diagnose(upper), config.planet);
  const auto lower_energy = mps::diagnose_dry_hydrostatic_budgets(
      driver.grid(), lower, driver.diagnose(lower), config.planet);
  const mps::Real measured =
      (upper_energy.total_energy_j - lower_energy.total_energy_j) / (2.0 * epsilon);
  const mps::Real expected = rhs.radiation_diagnostics.dry_thermal_energy_rate_w;
  MPS_CHECK_NEAR(measured, expected, 1e-6 * std::max(1.0, std::abs(expected)));
}

MPS_TEST_CASE("gray RHS evaluates orbit at the state time") {
  const auto config = gray_config();
  const mps::DryHydrostaticDriver driver(config);
  auto first = driver.initial_state();
  auto later = first;
  later.time_s = config.orbit->period_s / 4.0;
  const auto first_physics = driver.rhs_components(first).physics;
  const auto later_physics = driver.rhs_components(later).physics;
  bool differs = false;
  for (std::size_t cell = 0; cell < first.surface_temperature_k.size(); ++cell)
    differs = differs || first_physics.surface_temperature_k_s[cell] !=
                             later_physics.surface_temperature_k_s[cell];
  MPS_CHECK(differs);
}

MPS_TEST_CASE("explicit gray coupling conserves dry and tracer mass") {
  const auto config = gray_config();
  const mps::DryHydrostaticDriver driver(config);
  auto state = driver.initial_state();
  const mps::Real initial_mass = dry_mass(driver, state);
  const auto initial_tracer = state.tracer_mass_kg_m2;
  const auto initial_surface = state.surface_temperature_k;
  std::vector<mps::GrayRadiationBudget> budgets;
  driver.advance(
      state, config.run.end_time_s,
      [&](const mps::DryHydrostaticState& sampled, const mps::DryHydrostaticDerived*,
          const mps::DryHydrostaticStepDiagnostics& diagnostics) {
        if (sampled.step > 0) budgets.push_back(diagnostics.radiation_budget);
      });
  MPS_CHECK_NEAR(dry_mass(driver, state), initial_mass, 1e-13 * initial_mass);
  mps::Real initial_tracer_total = 0.0;
  mps::Real final_tracer_total = 0.0;
  for (std::size_t index = 0; index < initial_tracer.size(); ++index) {
    const std::size_t cell = index / static_cast<std::size_t>(config.vertical.levels);
    const mps::Real area = driver.grid().cells()[cell].area_m2;
    initial_tracer_total += area * initial_tracer[index];
    final_tracer_total += area * state.tracer_mass_kg_m2[index];
  }
  MPS_CHECK_NEAR(final_tracer_total, initial_tracer_total,
                 1e-13 * initial_tracer_total);
  MPS_CHECK(state.surface_temperature_k != initial_surface);
  mps::Real integrated_storage = 0.0;
  mps::Real integrated_interface_residual = 0.0;
  for (const auto& budget : budgets) {
    integrated_storage += budget.surface_storage_change_j;
    integrated_interface_residual += budget.interface_conservation_residual_j;
  }
  const mps::Real capacity = mps::mixed_surface_heat_capacity(
      config.surface->uniform_land_fraction, config.surface->land_heat_capacity_j_m2_k,
      config.surface->ocean_heat_capacity_j_m2_k);
  mps::Real measured_storage = 0.0;
  for (std::size_t cell = 0; cell < driver.grid().cell_count(); ++cell)
    measured_storage += driver.grid().cells()[cell].area_m2 * capacity *
                        (state.surface_temperature_k[cell] - initial_surface[cell]);
  MPS_CHECK_NEAR(integrated_storage, measured_storage,
                 1e-12 * std::max(1.0, std::abs(measured_storage)));
  MPS_CHECK(std::abs(integrated_interface_residual) <
            1e-12 * config.run.end_time_s *
                config.orbit->stellar_flux_at_semimajor_axis_w_m2 *
                driver.grid().total_area_m2());
}

MPS_TEST_CASE("semi implicit gray coupling advances with explicit source checks") {
  auto config = semi_implicit_config(10.0, 3);
  config.run.end_time_s = 10.0;
  config.validate();
  const mps::DryHydrostaticDriver driver(config);
  auto state = driver.initial_state();
  std::vector<mps::DryHydrostaticStepDiagnostics> diagnostics;
  driver.advance(state, config.run.end_time_s,
                 [&](const mps::DryHydrostaticState&, const mps::DryHydrostaticDerived*,
                     const mps::DryHydrostaticStepDiagnostics& sample) {
                   if (sample.accepted_time_step_s > 0.0) diagnostics.push_back(sample);
                 });
  MPS_CHECK_EQ(state.time_s, config.run.end_time_s);
  MPS_CHECK(!diagnostics.empty());
  MPS_CHECK(diagnostics.back().radiation_rates.toa_incoming_shortwave_power_w > 0.0);
}

MPS_TEST_CASE("centered gray semi implicit path has time and iteration convergence") {
  const auto coarse = run(semi_implicit_config(20.0, 5));
  const auto medium = run(semi_implicit_config(10.0, 5));
  const auto fine = run(semi_implicit_config(5.0, 5));
  const mps::Real coarse_difference = normalized_difference(coarse, medium, 2);
  const mps::Real fine_difference = normalized_difference(medium, fine, 2);
  MPS_CHECK(coarse_difference > 0.0);
  MPS_CHECK(fine_difference > 0.0);
  MPS_CHECK(std::log2(coarse_difference / fine_difference) >= 1.7);

  const auto three_iterations = run(semi_implicit_config(20.0, 3));
  const mps::Real iteration_error = normalized_difference(three_iterations, coarse, 2);
  MPS_CHECK(iteration_error < coarse_difference);
}

MPS_TEST_CASE("rejected radiation attempts do not enter accepted energy budgets") {
  auto config = gray_config();
  config.run.end_time_s = 0.1;
  config.run.time_step_s = 0.1;
  config.surface->land_heat_capacity_j_m2_k = 100.0;
  config.surface->ocean_heat_capacity_j_m2_k = 100.0;
  config.surface->air_exchange_coefficient_w_m2_k = 0.0;
  config.orbit->stellar_flux_at_semimajor_axis_w_m2 = 1e8;
  const mps::DryHydrostaticDriver driver(config);
  auto state = driver.initial_state();
  std::size_t retries = 0;
  mps::Real accepted_duration = 0.0;
  mps::Real incident_energy = 0.0;
  driver.advance(state, config.run.end_time_s,
                 [&](const mps::DryHydrostaticState&, const mps::DryHydrostaticDerived*,
                     const mps::DryHydrostaticStepDiagnostics& step) {
                   retries += step.retry_count;
                   MPS_CHECK(step.retry_count == step.cfl_retry_count +
                                                     step.invariant_retry_count +
                                                     step.solver_retry_count);
                   accepted_duration += step.accepted_time_step_s;
                   incident_energy +=
                       step.radiation_budget.toa_incoming_shortwave_energy_j;
                 });
  MPS_CHECK(retries > 0);
  MPS_CHECK_NEAR(accepted_duration, config.run.end_time_s, 1e-14);
  mps::DryHydrostaticRhs initial_rhs;
  driver.rhs(driver.initial_state(), initial_rhs);
  // Insolation changes negligibly over 0.1 seconds. Counting rejected attempts
  // would exceed this independent boundary-flux integral by orders of magnitude.
  const auto expected =
      config.run.end_time_s *
      initial_rhs.radiation_diagnostics.toa_incoming_shortwave_power_w;
  MPS_CHECK_NEAR(incident_energy, expected, 1e-8 * expected);
}

MPS_TEST_CASE("cancellation flushes an unsampled radiative interval once") {
  for (const bool implicit : {false, true}) {
    auto config = implicit ? semi_implicit_config(10.0, 5) : gray_config();
    config.diagnostics.interval_steps = 100;
    const mps::DryHydrostaticDriver driver(config);
    auto state = driver.initial_state();
    std::size_t accepted_callbacks = 0, samples = 0;
    double pending_energy = 0.0, written_energy = 0.0;
    driver.advance(
        state, config.run.end_time_s,
        [&](const mps::DryHydrostaticState&, const mps::DryHydrostaticDerived* derived,
            const mps::DryHydrostaticStepDiagnostics& step) {
          if (step.accepted_time_step_s > 0.0) ++accepted_callbacks;
          pending_energy += step.radiation_budget.toa_incoming_shortwave_energy_j;
          if (derived != nullptr) {
            ++samples;
            written_energy += pending_energy;
            pending_energy = 0.0;
          }
        },
        [&] { return state.step >= 1; });
    MPS_CHECK_EQ(state.step, 1U);
    MPS_CHECK_EQ(accepted_callbacks, 1U);
    MPS_CHECK_EQ(samples, 2U);
    MPS_CHECK_NEAR(pending_energy, 0.0, 0.0);
    MPS_CHECK(written_energy > 0.0);
  }
}

MPS_TEST_CASE("dry convection runs only after accepted explicit and implicit steps") {
  for (const bool implicit : {false, true}) {
    const auto config = convection_config(implicit);
    const mps::DryHydrostaticDriver driver(config);
    auto state = driver.initial_state();
    impose_unstable_columns(driver, state);
    const auto before_rhs = state.potential_temperature_mass_k_kg_m2;
    static_cast<void>(driver.rhs(state));
    MPS_CHECK(state.potential_temperature_mass_k_kg_m2 == before_rhs);

    mps::DryHydrostaticStepDiagnostics accepted;
    driver.advance(
        state, 10.0,
        [&](const mps::DryHydrostaticState&, const mps::DryHydrostaticDerived*,
            const mps::DryHydrostaticStepDiagnostics& step) {
          if (step.accepted_time_step_s > 0.0) accepted = step;
        });
    MPS_CHECK_EQ(state.step, 1U);
    MPS_CHECK(accepted.convection.adjusted_column_count > 0);
    MPS_CHECK(accepted.convection.unstable_interface_fraction_before > 0.0);
    MPS_CHECK_EQ(accepted.convection.unstable_interface_fraction_after, 0.0);
    MPS_CHECK(accepted.convection.minimum_theta_difference_after_k >= -1e-10);
    MPS_CHECK(accepted.convection_column_call_count >= driver.grid().cell_count());
    const auto adjusted = driver.diagnose(state);
    for (std::size_t cell = 0; cell < adjusted.cells; ++cell) {
      const auto upper = mps::dry_hydrostatic_offset(cell, 0, adjusted.levels);
      const auto lower = mps::dry_hydrostatic_offset(cell, 1, adjusted.levels);
      MPS_CHECK(adjusted.potential_temperature_k[upper] + 1e-10 >=
                adjusted.potential_temperature_k[lower]);
    }
  }
}

MPS_TEST_CASE("convective restart segmentation reproduces uninterrupted steps") {
  for (const bool implicit : {false, true}) {
    const auto config = convection_config(implicit);
    const mps::DryHydrostaticDriver full_driver(config);
    auto uninterrupted = full_driver.initial_state();
    impose_unstable_columns(full_driver, uninterrupted);
    auto segmented = uninterrupted;
    full_driver.advance(uninterrupted, config.run.end_time_s);
    full_driver.advance(segmented, 10.0);
    const mps::DryHydrostaticDriver restarted_driver(config);
    restarted_driver.advance(segmented, config.run.end_time_s);
    MPS_CHECK_NEAR(normalized_difference(uninterrupted, segmented, 2), 0.0, 1e-14);
  }
}

int main() { return mps::test::run_all(); }
