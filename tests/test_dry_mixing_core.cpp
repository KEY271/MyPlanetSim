#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "myplanetsim/dynamics/dry_hydrostatic_driver.hpp"
#include "support/test.hpp"

namespace {

[[nodiscard]] mps::ExperimentConfig mixing_config(const bool implicit) {
  mps::ExperimentConfig config{
      .kind = mps::ExperimentKind::kDryHydrostatic,
      .planet = {6.37122e6, 7.292115e-5, 9.80616, 287.0, 1004.0, 100000.0},
      .run = {0.0, 20.0, 10.0, 7},
      .grid = {2},
      .vertical = {.levels = 3,
                   .a_half_pa = {1000.0, 500.0, 250.0, 0.0},
                   .b_half = {0.0, 0.33, 0.66, 1.0},
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
      .convection = {.kind = mps::ConvectionKind::kDryAdjustment,
                     .stability_tolerance_k = 1e-10},
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
              .uniform_land_fraction = 0.25,
              .land_heat_capacity_j_m2_k = 2e6,
              .ocean_heat_capacity_j_m2_k = 2e8,
              .initial_temperature_k = 310.0,
              .albedo = 0.25,
              .emissivity = 0.95,
              .air_exchange_coefficient_w_m2_k = 0.0,
              .internal_heat_flux_w_m2 = 0.1,
              .land_roughness_momentum_m = 0.1,
              .land_roughness_heat_m = 0.01,
              .ocean_roughness_momentum_m = 0.001,
              .ocean_roughness_heat_m = 0.0001,
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
  if (implicit) {
    config.dry_hydrostatic.time_integrator =
        mps::DryHydrostaticTimeIntegrator::kSemiImplicit;
    config.dry_hydrostatic.advective_cfl = 0.45;
    config.semi_implicit = mps::SemiImplicitParameters{
        .reference_surface_pressure_pa = 100000.0,
        .reference_temperature_k = 280.0,
        .reference_update = mps::SemiImplicitReferenceUpdate::kFixed,
        .implicit_weight = 0.5,
        .wave_cfl_threshold = 0.45,
        .maximum_implicit_modes = 3,
        .nonlinear_iterations = 5,
        .linear_relative_tolerance = 1e-8,
        .linear_absolute_tolerance = 1e-12,
        .linear_maximum_iterations = 80,
        .gmres_restart = 20,
        .minimum_time_step_s = 0.01};
  }
  config.validate();
  return config;
}

[[nodiscard]] mps::Real normalized_difference(const mps::DryHydrostaticState& first,
                                              const mps::DryHydrostaticState& second,
                                              const std::size_t levels) {
  const auto first_values = mps::flatten_dry_hydrostatic_surface_state(first, levels);
  const auto second_values = mps::flatten_dry_hydrostatic_surface_state(second, levels);
  mps::Real result = 0.0;
  for (std::size_t index = 0; index < first_values.size(); ++index) {
    const mps::Real scale = std::max(1.0, std::abs(first_values[index]));
    result =
        std::max(result, std::abs(first_values[index] - second_values[index]) / scale);
  }
  return result;
}

}  // namespace

MPS_TEST_CASE("boundary layer replaces legacy gray surface exchange") {
  const auto config = mixing_config(false);
  const mps::DryHydrostaticDriver driver(config);
  const auto state = driver.initial_state();
  const auto rhs = driver.rhs(state);
  const auto components = driver.rhs_components(state);
  MPS_CHECK_EQ(rhs.radiation_diagnostics.sensible_to_atmosphere_power_w, 0.0);
  MPS_CHECK_EQ(rhs.radiation_diagnostics.rayleigh_drag_work_w, 0.0);
  for (const auto momentum : components.physics.tendency.momentum) {
    MPS_CHECK_EQ(momentum.x, 0.0);
    MPS_CHECK_EQ(momentum.y, 0.0);
    MPS_CHECK_EQ(momentum.z, 0.0);
  }
}

MPS_TEST_CASE("full dry mixing closes accepted-step budgets for both integrators") {
  for (const bool implicit : {false, true}) {
    const auto config = mixing_config(implicit);
    const mps::DryHydrostaticDriver driver(config);
    auto state = driver.initial_state();
    const auto initial_surface_temperature = state.surface_temperature_k;
    mps::DryHydrostaticStepDiagnostics accepted;
    driver.advance(
        state, 10.0,
        [&](const mps::DryHydrostaticState&, const mps::DryHydrostaticDerived*,
            const mps::DryHydrostaticStepDiagnostics& step) {
          if (step.accepted_time_step_s > 0.0) accepted = step;
        });
    MPS_CHECK_EQ(state.step, 1U);
    MPS_CHECK(accepted.boundary_layer_column_call_count >= driver.grid().cell_count());
    MPS_CHECK(accepted.convection_column_call_count >= driver.grid().cell_count());
    const mps::Real heat_scale =
        std::max(1.0, std::abs(accepted.boundary_layer.atmospheric_heat_change_j) +
                          std::abs(accepted.boundary_layer.surface_heat_change_j) +
                          accepted.boundary_layer.returned_dissipation_heat_j);
    MPS_CHECK_NEAR(accepted.boundary_layer.heat_budget_residual_j, 0.0,
                   1e-10 * heat_scale);
    MPS_CHECK(std::abs(accepted.boundary_layer.tracer_mass_change_kg) <= 1e-12);
    MPS_CHECK(accepted.boundary_layer.mean_boundary_layer_height_m > 0.0);
    MPS_CHECK(accepted.convection.unstable_interface_fraction_after == 0.0);
    MPS_CHECK_EQ(accepted.radiation_budget.sensible_to_atmosphere_energy_j, 0.0);
    MPS_CHECK_EQ(accepted.radiation_budget.rayleigh_drag_energy_j, 0.0);
    const mps::Real surface_capacity =
        mps::mixed_surface_heat_capacity(config.surface->uniform_land_fraction,
                                         config.surface->land_heat_capacity_j_m2_k,
                                         config.surface->ocean_heat_capacity_j_m2_k);
    mps::Real measured_surface_change = 0.0;
    for (std::size_t cell = 0; cell < driver.grid().cell_count(); ++cell)
      measured_surface_change +=
          driver.grid().cells()[cell].area_m2 * surface_capacity *
          (state.surface_temperature_k[cell] - initial_surface_temperature[cell]);
    const mps::Real attributed_surface_change =
        accepted.radiation_budget.surface_storage_change_j +
        accepted.boundary_layer.surface_heat_change_j;
    MPS_CHECK_NEAR(measured_surface_change, attributed_surface_change,
                   1e-8 * std::max(1.0, std::abs(measured_surface_change)));
  }
}

MPS_TEST_CASE("full mixing restart segmentation matches uninterrupted steps") {
  for (const bool implicit : {false, true}) {
    const auto config = mixing_config(implicit);
    const mps::DryHydrostaticDriver first_driver(config);
    auto uninterrupted = first_driver.initial_state();
    auto segmented = uninterrupted;
    first_driver.advance(uninterrupted, config.run.end_time_s);
    first_driver.advance(segmented, 10.0);
    const mps::DryHydrostaticDriver restarted_driver(config);
    restarted_driver.advance(segmented, config.run.end_time_s);
    MPS_CHECK_NEAR(normalized_difference(uninterrupted, segmented, 3), 0.0, 1e-14);
  }
}

MPS_TEST_CASE("post-dynamics mixing failure rolls back the complete candidate") {
  auto config = mixing_config(true);
  config.surface->land_roughness_momentum_m = 1e6;
  config.surface->land_roughness_heat_m = 1e6;
  config.surface->ocean_roughness_momentum_m = 1e6;
  config.surface->ocean_roughness_heat_m = 1e6;
  config.validate();
  const mps::DryHydrostaticDriver driver(config);
  auto state = driver.initial_state();
  const auto initial = state;
  MPS_CHECK_THROWS_AS(driver.advance(state, 10.0), std::invalid_argument);
  MPS_CHECK_NEAR(normalized_difference(state, initial, 3), 0.0, 0.0);
  MPS_CHECK_EQ(state.time_s, initial.time_s);
  MPS_CHECK_EQ(state.step, initial.step);
}

int main() { return mps::test::run_all(); }
