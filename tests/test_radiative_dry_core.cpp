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
  driver.advance(state, config.run.end_time_s);
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

int main() { return mps::test::run_all(); }
