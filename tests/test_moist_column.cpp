#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

#include "myplanetsim/physics/boundary_layer.hpp"
#include "myplanetsim/physics/dry_convective_adjustment.hpp"
#include "myplanetsim/physics/gray_radiation.hpp"
#include "myplanetsim/physics/saturation_adjustment.hpp"
#include "myplanetsim/physics/simple_betts_miller.hpp"
#include "myplanetsim/vertical/hybrid_pressure_coordinate.hpp"
#include "myplanetsim/vertical/hydrostatic_column.hpp"
#include "support/test.hpp"

namespace {

struct ColumnSummary {
  double mean_temperature_k = 0.0;
  double precipitable_water_kg_m2 = 0.0;
  double surface_temperature_k = 0.0;
  double water_residual_kg_m2 = 0.0;
  double energy_residual_j_m2 = 0.0;
  double maximum_local_heat_residual_j_m2 = 0.0;
  double maximum_supersaturation = 0.0;
};

[[nodiscard]] double moist_enthalpy(const std::vector<double>& temperature,
                                    const std::vector<double>& vapor,
                                    const std::vector<double>& mass,
                                    const double surface_temperature) {
  constexpr mps::DiluteMoistThermodynamics moist;
  double result = 4.0e7 * surface_temperature;
  for (std::size_t level = 0; level < mass.size(); ++level)
    result += mass[level] * (moist.heat_capacity_cp_j_kg_k * temperature[level] +
                             moist.latent_heat_vaporization_j_kg * vapor[level]);
  return result;
}

[[nodiscard]] ColumnSummary run_column(const std::size_t levels, const double dt,
                                       const double final_time = 3600.0) {
  constexpr double gravity = 9.80616;
  constexpr double gas_constant = 287.0;
  constexpr double cp = 1004.0;
  constexpr double reference_pressure = 100000.0;
  constexpr mps::DiluteMoistThermodynamics moist;
  const auto coefficients =
      mps::surface_refined_sigma_coefficients(1000.0, levels, 5.0);
  const mps::AtmosphericHybridCoordinate coordinate(coefficients, 90000.0, 110000.0,
                                                    10.0);
  const auto geometry =
      coordinate.geometry(100000.0, gravity, gas_constant, cp, reference_pressure);
  const auto& mass = geometry.air_mass_kg_m2;
  std::vector<double> temperature(levels);
  std::vector<double> theta(levels);
  std::vector<double> vapor(levels);
  std::vector<mps::Vec3> velocity(levels, {3.0, 0.0, 0.0});
  for (std::size_t level = 0; level < levels; ++level) {
    const double sigma = geometry.pressure_full_pa[level] / reference_pressure;
    temperature[level] = 220.0 + 70.0 * std::pow(sigma, 0.25);
    theta[level] = temperature[level] / geometry.exner_full[level];
    vapor[level] =
        0.55 * mps::saturation_mixing_ratio(temperature[level],
                                            geometry.pressure_full_pa[level], moist);
  }
  double surface_temperature = 294.0;
  double ocean_water_ledger = 0.0;
  double precipitation_ledger = 0.0;
  double integrated_external_heat = 0.0;
  const double initial_water = [&] {
    double value = 0.0;
    for (std::size_t level = 0; level < levels; ++level)
      value += mass[level] * vapor[level];
    return value;
  }();
  const double initial_enthalpy =
      moist_enthalpy(temperature, vapor, mass, surface_temperature);
  mps::BoundaryLayerColumnWorkspace boundary_workspace;
  mps::BoundaryLayerColumnResult boundary_result;
  mps::DryConvectiveAdjustmentWorkspace dry_workspace;
  mps::DryConvectiveAdjustmentResult dry_result;
  mps::SimpleBettsMillerResult moist_result;
  mps::GrayRadiativeColumnWorkspace radiation_workspace;
  mps::GrayRadiativeColumnTendency radiation;
  double maximum_water_residual = 0.0;
  double maximum_process_heat_residual = 0.0;
  double maximum_supersaturation = 0.0;
  const std::size_t steps = static_cast<std::size_t>(std::llround(final_time / dt));
  for (std::size_t step = 0; step < steps; ++step) {
    mps::gray_radiative_column_tendency(
        {.radiation = {.pressure_half_pa = geometry.pressure_half_pa,
                       .temperature_k = temperature,
                       .surface_temperature_k = surface_temperature,
                       .gravity_m_s2 = gravity,
                       .stellar_flux_w_m2 = 1000.0,
                       .cosine_solar_zenith = 0.5,
                       .surface_albedo = 0.3,
                       .surface_emissivity = 1.0,
                       .parameters = {.shortwave_absorption_m2_kg = 0.0,
                                      .longwave_absorption_ref_m2_kg = 3e-4,
                                      .reference_pressure_pa = reference_pressure,
                                      .longwave_pressure_exponent = 1.0,
                                      .longwave_diffusivity_factor = 1.66,
                                      .shortwave_diffuse_factor = 1.66,
                                      .cfl = 0.5}},
         .exner_full = geometry.exner_full,
         .heat_capacity_cp_j_kg_k = cp,
         .surface_heat_capacity_j_m2_k = 4.0e7,
         .air_exchange_coefficient_w_m2_k = 0.0,
         .internal_heat_flux_w_m2 = 0.0},
        radiation, radiation_workspace);
    integrated_external_heat -= dt * radiation.fluxes.net_flux_w_m2.front();
    for (std::size_t level = 0; level < levels; ++level) {
      theta[level] +=
          dt * radiation.potential_temperature_mass_k_kg_m2_s[level] / mass[level];
      temperature[level] = theta[level] * geometry.exner_full[level];
    }
    surface_temperature += dt * radiation.surface_temperature_k_s;

    const auto hydrostatic =
        mps::integrate_hydrostatic_column(geometry, theta, cp, gravity, 0.0);
    std::vector<double> height_half = hydrostatic.height_half_m;
    std::vector<double> height_full = hydrostatic.height_full_m;
    height_half.back() = 0.0;
    const auto bulk = mps::diagnose_boundary_layer_column(
        {.potential_temperature_k = theta,
         .temperature_k = temperature,
         .velocity_m_s = velocity,
         .pressure_half_pa = geometry.pressure_half_pa,
         .exner_half = geometry.exner_half,
         .height_half_m = height_half,
         .height_full_m = height_full,
         .surface_temperature_k = surface_temperature,
         .surface_exner = geometry.exner_half.back(),
         .gravity_m_s2 = gravity,
         .gas_constant_j_kg_k = gas_constant,
         .heat_capacity_cp_j_kg_k = cp,
         .critical_richardson = 1.0,
         .turbulent_prandtl = 1.0,
         .gustiness_m_s = 1.0,
         .land_fraction = 0.0,
         .land_roughness = {.momentum_m = 0.1, .heat_m = 0.01},
         .ocean_roughness = {.momentum_m = 0.001, .heat_m = 0.0001}});
    mps::implicit_boundary_layer_column(
        {.potential_temperature_k = theta,
         .velocity_m_s = velocity,
         .tracer_mixing_ratio = vapor,
         .air_mass_kg_m2 = mass,
         .exner_full = geometry.exner_full,
         .exner_half = geometry.exner_half,
         .height_full_m = height_full,
         .density_half_kg_m3 = bulk.density_half_kg_m3,
         .eddy_diffusivity_momentum_m2_s = bulk.eddy_diffusivity_momentum_m2_s,
         .eddy_diffusivity_heat_m2_s = bulk.eddy_diffusivity_heat_m2_s,
         .eddy_diffusivity_tracer_m2_s = bulk.eddy_diffusivity_tracer_m2_s,
         .surface_temperature_k = surface_temperature,
         .surface_exner = geometry.exner_half.back(),
         .surface_heat_capacity_j_m2_k = 4.0e7,
         .surface_heat_conductance_w_m2_k = bulk.surface_heat_conductance_w_m2_k,
         .surface_drag_conductance_kg_m2_s = bulk.surface_drag_conductance_kg_m2_s,
         .heat_capacity_cp_j_kg_k = cp,
         .time_step_s = dt,
         .enable_surface_water_exchange = true,
         .surface_pressure_pa = 100000.0,
         .land_fraction = 0.0,
         .land_water_kg_m2 = 0.0,
         .bucket_capacity_kg_m2 = 0.0,
         .surface_water_conductance_land_kg_m2_s = 0.0,
         .surface_water_conductance_ocean_kg_m2_s =
             bulk.surface_water_conductance_ocean_kg_m2_s,
         .moist_thermodynamics = moist},
        boundary_result, boundary_workspace);
    theta = boundary_result.potential_temperature_k;
    velocity = boundary_result.velocity_m_s;
    vapor = boundary_result.tracer_mixing_ratio;
    surface_temperature = boundary_result.surface_temperature_k;
    ocean_water_ledger += boundary_result.diagnostics.ocean_water_change_kg_m2;
    maximum_water_residual =
        std::max(maximum_water_residual,
                 std::abs(boundary_result.diagnostics.water_budget_residual_kg_m2));
    maximum_process_heat_residual = std::max(
        maximum_process_heat_residual,
        std::abs(boundary_result.diagnostics.moist_enthalpy_budget_residual_j_m2));

    mps::dry_convective_adjustment({.potential_temperature_k = theta,
                                    .air_mass_kg_m2 = mass,
                                    .exner_full = geometry.exner_full,
                                    .exner_half = geometry.exner_half,
                                    .heat_capacity_cp_j_kg_k = cp,
                                    .heat_capacity_cv_j_kg_k = cp - gas_constant,
                                    .stability_tolerance_k = 1e-10},
                                   dry_result, dry_workspace);
    theta = dry_result.adjusted_potential_temperature_k;
    for (std::size_t level = 0; level < levels; ++level)
      temperature[level] = theta[level] * geometry.exner_full[level];
    mps::simple_betts_miller_adjustment({.temperature_k = temperature,
                                         .vapor_mixing_ratio = vapor,
                                         .air_mass_kg_m2 = mass,
                                         .pressure_full_pa = geometry.pressure_full_pa,
                                         .pressure_half_pa = geometry.pressure_half_pa,
                                         .gravity_m_s2 = gravity,
                                         .relative_humidity_reference = 0.8,
                                         .relaxation_time_s = 7200.0,
                                         .time_step_s = dt,
                                         .minimum_temperature_k = 150.0,
                                         .thermodynamics = moist},
                                        moist_result);
    temperature = moist_result.temperature_k;
    vapor = moist_result.vapor_mixing_ratio;
    precipitation_ledger += moist_result.diagnostics.convective_rain_kg_m2;
    maximum_process_heat_residual =
        std::max(maximum_process_heat_residual,
                 std::abs(moist_result.diagnostics.moist_enthalpy_change_j_m2));
    std::vector<double> theta_mass(levels);
    std::vector<double> vapor_mass(levels);
    for (std::size_t level = 0; level < levels; ++level) {
      theta[level] = temperature[level] / geometry.exner_full[level];
      theta_mass[level] = mass[level] * theta[level];
      vapor_mass[level] = mass[level] * vapor[level];
    }
    const auto saturation = mps::adjust_saturation_column(
        mass, geometry.pressure_full_pa, geometry.exner_full, theta_mass, vapor_mass,
        moist);
    precipitation_ledger += saturation.condensed_water_kg_m2;
    maximum_process_heat_residual = std::max(maximum_process_heat_residual,
                                             std::abs(saturation.enthalpy_change_j_m2));
    maximum_supersaturation =
        std::max(maximum_supersaturation, saturation.maximum_supersaturation);
    for (std::size_t level = 0; level < levels; ++level) {
      theta[level] = theta_mass[level] / mass[level];
      temperature[level] = theta[level] * geometry.exner_full[level];
      vapor[level] = vapor_mass[level] / mass[level];
    }
  }
  double final_water = ocean_water_ledger + precipitation_ledger;
  double mean_temperature = 0.0;
  double total_air_mass = 0.0;
  for (std::size_t level = 0; level < levels; ++level) {
    final_water += mass[level] * vapor[level];
    mean_temperature += mass[level] * temperature[level];
    total_air_mass += mass[level];
  }
  const double final_enthalpy =
      moist_enthalpy(temperature, vapor, mass, surface_temperature);
  return {.mean_temperature_k = mean_temperature / total_air_mass,
          .precipitable_water_kg_m2 =
              final_water - ocean_water_ledger - precipitation_ledger,
          .surface_temperature_k = surface_temperature,
          .water_residual_kg_m2 =
              std::max(maximum_water_residual, std::abs(final_water - initial_water)),
          .energy_residual_j_m2 =
              std::abs(final_enthalpy - initial_enthalpy - integrated_external_heat),
          .maximum_local_heat_residual_j_m2 = maximum_process_heat_residual,
          .maximum_supersaturation = maximum_supersaturation};
}

[[nodiscard]] double distance(const ColumnSummary& first, const ColumnSummary& second) {
  return std::hypot(first.mean_temperature_k - second.mean_temperature_k,
                    first.surface_temperature_k - second.surface_temperature_k) +
         std::abs(first.precipitable_water_kg_m2 - second.precipitable_water_kg_m2);
}

}  // namespace

MPS_TEST_CASE("independent radiative moist column closes water and local heat") {
  const auto result = run_column(20, 75.0);
  MPS_CHECK(result.water_residual_kg_m2 <= 1e-10);
  MPS_CHECK_NEAR(result.maximum_local_heat_residual_j_m2, 0.0, 1e-3);
  // The composed radiation/physics ledger closes to much better than 1e-8 of
  // the O(1e10 J m-2) column reservoir; each individual moist process is checked
  // against the tighter absolute tolerance above.
  MPS_CHECK(result.energy_residual_j_m2 <= 200.0);
  MPS_CHECK(result.maximum_supersaturation <= 1e-12);
  MPS_CHECK(result.precipitable_water_kg_m2 > 0.0);
}

MPS_TEST_CASE("smooth moist column coupling is first order in time") {
  const auto reference = run_column(20, 18.75);
  const std::array<double, 3> errors{distance(run_column(20, 300.0), reference),
                                     distance(run_column(20, 150.0), reference),
                                     distance(run_column(20, 75.0), reference)};
  MPS_CHECK(std::log2(errors[0] / errors[1]) >= 0.8);
  MPS_CHECK(std::log2(errors[1] / errors[2]) >= 0.8);
}

MPS_TEST_CASE("BL candidate intervals retain moist column constraints") {
  for (const double dt : {75.0, 150.0, 300.0, 600.0, 900.0}) {
    const auto result = run_column(20, dt);
    MPS_CHECK(std::isfinite(result.mean_temperature_k));
    MPS_CHECK(std::isfinite(result.surface_temperature_k));
    MPS_CHECK(result.precipitable_water_kg_m2 > 0.0);
    MPS_CHECK(result.water_residual_kg_m2 <= 1e-10);
    MPS_CHECK(result.maximum_supersaturation <= 1e-12);
    MPS_CHECK_NEAR(result.maximum_local_heat_residual_j_m2, 0.0, 1e-3);
  }
}

MPS_TEST_CASE("surface-refined moist columns remain bounded at K20 K40 K80") {
  const auto k20 = run_column(20, 75.0);
  const auto k40 = run_column(40, 75.0);
  const auto k80 = run_column(80, 75.0);
  for (const auto& result : {k20, k40, k80}) {
    MPS_CHECK(std::isfinite(result.mean_temperature_k));
    MPS_CHECK(std::isfinite(result.precipitable_water_kg_m2));
    MPS_CHECK(result.precipitable_water_kg_m2 > 0.0);
    MPS_CHECK(result.water_residual_kg_m2 <= 1e-10);
  }
  MPS_CHECK(distance(k40, k80) < distance(k20, k40));
}

int main() { return mps::test::run_all(); }
