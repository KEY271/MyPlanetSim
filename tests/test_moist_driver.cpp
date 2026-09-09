#include <cmath>
#include <sstream>

#include "myplanetsim/config/experiment_config.hpp"
#include "myplanetsim/dynamics/dry_hydrostatic_driver.hpp"
#include "support/test.hpp"

namespace {

[[nodiscard]] mps::ExperimentConfig moist_config(const bool semi_implicit) {
  mps::ExperimentConfig config{
      .kind = mps::ExperimentKind::kDryHydrostatic,
      .planet = {6.37122e6, 7.292115e-5, 9.80616, 287.0, 1004.0, 100000.0},
      .run = {0.0, 600.0, 600.0, 7},
      .grid = {2},
      .vertical = {.levels = 2,
                   .a_half_pa = {20000.0, 10000.0, 0.0},
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
      .tracers = {{.name = "passive",
                   .role = mps::TracerRole::kPassive,
                   .initial_mixing_ratio = 0.1,
                   .require_nonnegative = true,
                   .horizontal_diffusion = true},
                  {.name = "water_vapor",
                   .role = mps::TracerRole::kWaterVapor,
                   .initial_relative_humidity = 0.6,
                   .require_nonnegative = true,
                   .horizontal_diffusion = false}},
      .physics = {.kind = mps::PhysicsKind::kGrayRadiation},
      .convection = {.kind = mps::ConvectionKind::kSimpleBettsMiller,
                     .stability_tolerance_k = 1e-10,
                     .relaxation_time_s = 7200.0,
                     .reference_relative_humidity = 0.8},
      .moisture = {.kind = mps::MoistureKind::kDiluteWater,
                   .condensation = mps::CondensationKind::kSaturationAdjustment,
                   .surface_exchange = mps::SurfaceMoistureExchange::kBulk,
                   .maximum_physics_substep_s = 300.0},
      .boundary_layer = {.kind = mps::BoundaryLayerKind::kBulkKProfile,
                         .integrator = mps::BoundaryLayerIntegrator::kBackwardEuler,
                         .critical_richardson = 1.0,
                         .turbulent_prandtl = 1.0,
                         .gustiness_m_s = 1.0},
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
              .uniform_land_fraction = 0.5,
              .land_heat_capacity_j_m2_k = 2e6,
              .ocean_heat_capacity_j_m2_k = 4e7,
              .initial_temperature_k = 310.0,
              .albedo = 0.3,
              .emissivity = 1.0,
              .air_exchange_coefficient_w_m2_k = 0.0,
              .internal_heat_flux_w_m2 = 0.0,
              .land_roughness_momentum_m = 0.1,
              .land_roughness_heat_m = 0.01,
              .ocean_roughness_momentum_m = 0.001,
              .ocean_roughness_heat_m = 0.0001,
              .hydrology_kind = mps::SurfaceHydrologyKind::kBucket,
              .hydrology_capacity_kg_m2 = 150.0,
              .hydrology_initial_fraction = 0.5,
              .cfl = 0.5,
          },
      .radiation =
          mps::RadiationParameters{
              .shortwave_absorption_m2_kg = 0.0,
              .longwave_absorption_ref_m2_kg = 3e-4,
              .reference_pressure_pa = 100000.0,
              .longwave_pressure_exponent = 1.0,
              .longwave_diffusivity_factor = 1.66,
              .shortwave_diffuse_factor = 1.66,
              .cfl = 0.5,
          },
  };
  if (semi_implicit) {
    config.dry_hydrostatic.time_integrator =
        mps::DryHydrostaticTimeIntegrator::kSemiImplicit;
    config.dry_hydrostatic.advective_cfl = 0.45;
    config.semi_implicit = mps::SemiImplicitParameters{
        .reference_surface_pressure_pa = 100000.0,
        .reference_temperature_k = 280.0,
        .reference_update = mps::SemiImplicitReferenceUpdate::kFixed,
        .implicit_weight = 0.5,
        .wave_cfl_threshold = 0.45,
        .maximum_implicit_modes = 2,
        .nonlinear_iterations = 3,
        .linear_relative_tolerance = 1e-8,
        .linear_absolute_tolerance = 1e-12,
        .linear_maximum_iterations = 80,
        .gmres_restart = 20,
        .minimum_time_step_s = 0.01,
    };
  }
  config.validate();
  return config;
}

[[nodiscard]] double water_inventory(const mps::DryHydrostaticDriver& driver,
                                     const mps::DryHydrostaticState& state) {
  const auto derived = driver.diagnose(state);
  const std::size_t volume = derived.cells * derived.levels;
  double total =
      state.cumulative_ocean_water_change_kg + state.cumulative_external_outflow_kg;
  for (std::size_t cell = 0; cell < derived.cells; ++cell) {
    const double area = driver.grid().cells()[cell].area_m2;
    for (std::size_t level = 0; level < derived.levels; ++level)
      total +=
          area * state.tracer_mass_kg_m2[volume + mps::dry_hydrostatic_offset(
                                                      cell, level, derived.levels)];
    total += area * driver.surface_boundary()->land_fraction()[cell] *
             state.land_water_kg_m2[cell];
  }
  return total;
}

}  // namespace

MPS_TEST_CASE("moist config round trips without changing its schema") {
  const auto config = moist_config(false);
  std::ostringstream encoded;
  mps::write_experiment_config(encoded, config);
  std::istringstream input(encoded.str());
  const auto decoded = mps::parse_experiment_config(input);
  MPS_CHECK(decoded.moisture.kind == mps::MoistureKind::kDiluteWater);
  MPS_CHECK(decoded.convection.kind == mps::ConvectionKind::kSimpleBettsMiller);
  MPS_CHECK(decoded.tracers[1].initial_relative_humidity.has_value());
  std::ostringstream second;
  mps::write_experiment_config(second, decoded);
  MPS_CHECK_EQ(second.str(), encoded.str());
}

MPS_TEST_CASE("moist global runs require gray radiation and a surface") {
  auto config = moist_config(false);
  config.physics.kind = mps::PhysicsKind::kNone;
  config.surface.reset();
  config.radiation.reset();
  MPS_CHECK_THROWS_AS(config.validate(), std::invalid_argument);
}

MPS_TEST_CASE("explicit and semi-implicit drivers accept subcycled moist steps") {
  for (const bool semi_implicit : {false, true}) {
    const auto config = moist_config(semi_implicit);
    const mps::DryHydrostaticDriver driver(config);
    auto state = driver.initial_state();
    const double initial_water = water_inventory(driver, state);
    std::size_t physics_substeps = 0;
    double diagnosed_evaporation = 0.0;
    driver.advance(state, config.run.end_time_s,
                   [&](const auto&, const auto*, const auto& step) {
                     physics_substeps += step.physics_substep_count;
                     diagnosed_evaporation += step.moisture.evaporation_kg;
                   });
    const double final_water = water_inventory(driver, state);
    MPS_CHECK_EQ(state.time_s, config.run.end_time_s);
    MPS_CHECK(physics_substeps >= 2U);
    MPS_CHECK(state.cumulative_evaporation_kg > 0.0);
    MPS_CHECK_NEAR(diagnosed_evaporation, state.cumulative_evaporation_kg,
                   1e-14 * state.cumulative_evaporation_kg);
    MPS_CHECK_NEAR(final_water, initial_water, 2e-12 * initial_water);
    for (const double water : state.land_water_kg_m2)
      MPS_CHECK(water >= 0.0 && water <= 150.0);
  }
}

MPS_TEST_CASE("process scheduler applies independent accepted intervals") {
  for (const bool semi_implicit : {false, true}) {
    auto config = moist_config(semi_implicit);
    config.physics_schedule = {
        .kind = mps::PhysicsScheduleKind::kProcessIntervals,
        .boundary_layer_maximum_update_interval_s = 200.0,
        .convection_diagnostic_interval_s = 300.0,
        .convection_update_mode = mps::ConvectionUpdateMode::kIntermittent,
        .radiation_diagnostic_interval_s = 600.0};
    config.validate();
    const mps::DryHydrostaticDriver driver(config);
    auto state = driver.initial_state();
    const double initial_water = water_inventory(driver, state);
    mps::DryHydrostaticStepDiagnostics accepted;
    driver.advance(state, config.run.end_time_s,
                   [&](const auto& sampled, const auto*, const auto& step) {
                     if (sampled.step > 0) accepted = step;
                   });
    MPS_CHECK_EQ(accepted.physics_substep_count, 4U);
    MPS_CHECK_EQ(accepted.radiation_column_call_count, driver.grid().cell_count());
    MPS_CHECK_EQ(accepted.boundary_layer_column_call_count,
                 3U * driver.grid().cell_count());
    MPS_CHECK_EQ(accepted.convection_column_call_count,
                 4U * driver.grid().cell_count());
    MPS_CHECK_EQ(accepted.physics_retry_count, 0U);
    MPS_CHECK(std::isfinite(accepted.radiation_budget.toa_net_upward_energy_j));
    MPS_CHECK_NEAR(water_inventory(driver, state), initial_water,
                   2e-12 * initial_water);
  }
}

MPS_TEST_CASE("moist checkpoint preserves bucket and accepted ledgers") {
  const auto config = moist_config(false);
  const mps::DryHydrostaticDriver driver(config);
  auto state = driver.initial_state();
  driver.advance(state, config.run.end_time_s);
  const auto flat = mps::flatten_dry_hydrostatic_moist_state(state, 2);
  const auto restored = mps::unflatten_dry_hydrostatic_moist_state(
      state.time_s, state.step, flat, driver.grid().cell_count(), 2, 2);
  MPS_CHECK(restored.land_water_kg_m2 == state.land_water_kg_m2);
  MPS_CHECK_EQ(restored.cumulative_evaporation_kg, state.cumulative_evaporation_kg);
  MPS_CHECK_EQ(restored.cumulative_ocean_water_change_kg,
               state.cumulative_ocean_water_change_kg);
  MPS_CHECK(restored.tracer_mass_kg_m2 == state.tracer_mass_kg_m2);
}

MPS_TEST_CASE("moist restart follows the uninterrupted accepted-step partition") {
  for (const bool semi_implicit : {false, true}) {
    for (const bool process_intervals : {false, true}) {
      auto config = moist_config(semi_implicit);
      config.run.end_time_s = 1200.0;
      if (process_intervals)
        config.physics_schedule = {
            .kind = mps::PhysicsScheduleKind::kProcessIntervals,
            .boundary_layer_maximum_update_interval_s = 300.0,
            .convection_diagnostic_interval_s = 600.0,
            .convection_update_mode = mps::ConvectionUpdateMode::kIntermittent,
            .radiation_diagnostic_interval_s = 600.0};
      config.validate();

      const mps::DryHydrostaticDriver continuous_driver(config);
      auto continuous = continuous_driver.initial_state();
      continuous_driver.advance(continuous, 1200.0);

      const mps::DryHydrostaticDriver first_driver(config);
      auto split = first_driver.initial_state();
      first_driver.advance(split, 600.0);
      const auto checkpoint = mps::flatten_dry_hydrostatic_moist_state(split, 2);
      auto restarted = mps::unflatten_dry_hydrostatic_moist_state(
          split.time_s, split.step, checkpoint, first_driver.grid().cell_count(), 2, 2);
      const mps::DryHydrostaticDriver restart_driver(config);
      restart_driver.advance(restarted, 1200.0);

      MPS_CHECK_EQ(restarted.time_s, continuous.time_s);
      MPS_CHECK_EQ(restarted.step, continuous.step);
      MPS_CHECK(mps::flatten_dry_hydrostatic_moist_state(restarted, 2) ==
                mps::flatten_dry_hydrostatic_moist_state(continuous, 2));
    }
  }
}

MPS_TEST_CASE("initial moisture outside the dilute domain is rejected") {
  auto config = moist_config(false);
  config.vertical.a_half_pa = {1000.0, 500.0, 0.0};
  config.validate();
  const mps::DryHydrostaticDriver driver(config);
  MPS_CHECK_THROWS_AS(driver.initial_state(), std::domain_error);
}

int main() { return mps::test::run_all(); }
