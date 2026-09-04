#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <numbers>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "myplanetsim/diagnostics/vertical_column_diagnostics.hpp"
#include "myplanetsim/dynamics/vertical_column_driver.hpp"
#include "myplanetsim/io/checkpoint.hpp"
#include "myplanetsim/io/run_metadata.hpp"
#include "myplanetsim/thermodynamics/dry_thermodynamics.hpp"
#include "myplanetsim/vertical/hydrostatic_column.hpp"
#include "myplanetsim/vertical/vertical_transport.hpp"
#include "support/test.hpp"

namespace {

mps::AtmosphericHybridCoordinate coordinate() {
  return {{{1000.0, 500.0, 0.0}, {0.0, 0.5, 1.0}}, 90000.0, 110000.0, 100.0};
}

mps::ExperimentConfig column_config() {
  mps::ExperimentConfig config{};
  config.kind = mps::ExperimentKind::kVerticalColumn;
  config.planet = mps::PlanetParameters::earth_like();
  config.run = {
      .start_time_s = 0.0, .end_time_s = 20.0, .time_step_s = 0.25, .random_seed = 0};
  config.vertical = {.test_case = mps::VerticalTestCase::kMovingSurfacePressure,
                     .levels = 2,
                     .a_half_pa = {1000.0, 500.0, 0.0},
                     .b_half = {0.0, 0.5, 1.0},
                     .surface_pressure_pa = 100000.0,
                     .minimum_surface_pressure_pa = 90000.0,
                     .maximum_surface_pressure_pa = 110000.0,
                     .minimum_pressure_thickness_pa = 100.0,
                     .surface_geopotential_m2_s2 = 0.0,
                     .initial_temperature_k = 280.0,
                     .initial_potential_temperature_k = 300.0,
                     .temperature_floor_k = 100.0,
                     .transport_scheme = mps::VerticalTransportScheme::kDonorCell,
                     .limiter = mps::VerticalLimiterKind::kNone,
                     .cfl = 0.5,
                     .forcing_amplitude = 5.0};
  config.diagnostics.interval_steps = 1;
  config.output_directory = "output";
  return config;
}

double spatial_transport_error(const std::size_t levels,
                               const mps::VerticalTransportScheme scheme) {
  const double layer_mass = 1.0 / static_cast<double>(levels);
  std::vector<mps::Real> air_mass(levels, layer_mass);
  std::vector<mps::Real> scalar_mass(levels);
  std::vector<mps::Real> horizontal_tendency(levels, 0.0);
  for (std::size_t k = 0; k < levels; ++k) {
    const double center = (static_cast<double>(k) + 0.5) * layer_mass;
    const double scalar = 2.0 + std::sin(2.0 * std::numbers::pi * center);
    scalar_mass[k] = layer_mass * scalar;
  }
  mps::VerticalMassFlux flux;
  flux.interface_flux_kg_m2_s.assign(levels + 1, 0.0);
  for (std::size_t interface = 1; interface < levels; ++interface) {
    const double x = static_cast<double>(interface) * layer_mass;
    flux.interface_flux_kg_m2_s[interface] = (x - 0.5) * std::sin(std::numbers::pi * x);
  }
  const auto numerical =
      mps::vertical_scalar_rhs(scalar_mass, air_mass, horizontal_tendency, flux, scheme,
                               mps::VerticalLimiterKind::kNone);
  double squared_error = 0.0;
  for (std::size_t k = 0; k < levels; ++k) {
    const double upper = static_cast<double>(k) * layer_mass;
    const double lower = static_cast<double>(k + 1) * layer_mass;
    const auto exact_flux = [](const double x) {
      if (x == 0.0 || x == 1.0) return 0.0;
      return (x - 0.5) * std::sin(std::numbers::pi * x) *
             (2.0 + std::sin(2.0 * std::numbers::pi * x));
    };
    const double exact = exact_flux(upper) - exact_flux(lower);
    const double error = (numerical[k] - exact) / layer_mass;
    squared_error += layer_mass * error * error;
  }
  return std::sqrt(squared_error);
}

double surface_pressure_quarter_period_error(const double time_step_s) {
  auto config = column_config();
  config.run = {.start_time_s = 10.0,
                .end_time_s = 26.0,
                .time_step_s = time_step_s,
                .random_seed = 0};
  config.vertical.forcing_amplitude = 100.0;
  config.vertical.minimum_surface_pressure_pa = 99000.0;
  config.vertical.maximum_surface_pressure_pa = 101000.0;
  const auto steps = static_cast<std::uint64_t>(4.0 / time_step_s);
  const auto result = mps::run_vertical_column(config, std::nullopt, steps);
  const double exact =
      config.vertical.surface_pressure_pa +
      config.vertical.forcing_amplitude * 16.0 / (2.0 * std::numbers::pi);
  return std::abs(result.state.surface_pressure_pa - exact);
}

}  // namespace

MPS_TEST_CASE("hybrid geometry telescopes mass and has interior full levels") {
  const auto geometry = coordinate().geometry(100000.0, 10.0, 287.0, 1004.0, 100000.0);
  MPS_CHECK_NEAR(std::accumulate(geometry.air_mass_kg_m2.begin(),
                                 geometry.air_mass_kg_m2.end(), 0.0),
                 9900.0, 1.0e-12);
  MPS_CHECK(geometry.pressure_full_pa[0] > geometry.pressure_half_pa[0]);
  MPS_CHECK(geometry.pressure_full_pa[0] < geometry.pressure_half_pa[1]);
  MPS_CHECK_THROWS_AS(mps::AtmosphericHybridCoordinate({{1000.0, 0.0}, {0.0, 0.9}},
                                                       90000.0, 110000.0, 100.0),
                      std::invalid_argument);
}

MPS_TEST_CASE("hybrid coordinate rejects malformed and unsafe profiles") {
  MPS_CHECK_THROWS_AS(mps::AtmosphericHybridCoordinate({{1000.0, 0.0}, {0.0}}, 90000.0,
                                                       110000.0, 100.0),
                      std::invalid_argument);
  MPS_CHECK_THROWS_AS(mps::AtmosphericHybridCoordinate(
                          {{std::numeric_limits<double>::quiet_NaN(), 0.0}, {0.0, 1.0}},
                          90000.0, 110000.0, 100.0),
                      std::invalid_argument);
  MPS_CHECK_THROWS_AS(mps::AtmosphericHybridCoordinate({{0.0, 0.0}, {0.0, 1.0}},
                                                       90000.0, 110000.0, 100.0),
                      std::invalid_argument);
  MPS_CHECK_THROWS_AS(
      mps::AtmosphericHybridCoordinate({{1000.0, 999.0, 0.0}, {0.0, 0.0, 1.0}}, 90000.0,
                                       110000.0, 100.0),
      std::invalid_argument);
  MPS_CHECK_THROWS_AS(
      mps::AtmosphericHybridCoordinate({{1000.0, 500.0, 0.0}, {0.0, 0.75, 0.5}},
                                       90000.0, 110000.0, 100.0),
      std::invalid_argument);
  for (const double surface_pressure : {90000.0, 110000.0}) {
    const auto geometry =
        coordinate().geometry(surface_pressure, 10.0, 287.0, 1004.0, 100000.0);
    MPS_CHECK(geometry.pressure_full_pa.front() > geometry.pressure_half_pa.front());
    MPS_CHECK(geometry.pressure_full_pa.back() < geometry.pressure_half_pa.back());
  }
}

MPS_TEST_CASE(
    "dry thermodynamics round trips and hydrostatic column is exact per layer") {
  const auto geometry = coordinate().geometry(100000.0, 10.0, 287.0, 1004.0, 100000.0);
  const auto theta = mps::thermodynamics::potential_temperature_from_temperature(
      280.0, geometry.pressure_full_pa[0], 287.0, 1004.0, 100000.0);
  MPS_CHECK_NEAR(mps::thermodynamics::temperature_from_potential_temperature(
                     theta, geometry.pressure_full_pa[0], 287.0, 1004.0, 100000.0),
                 280.0, 1.0e-12);
  const std::vector<mps::Real> values{300.0, 300.0};
  const auto hydrostatic =
      mps::integrate_hydrostatic_column(geometry, values, 1004.0, 10.0, 123.0);
  MPS_CHECK_NEAR(hydrostatic.geopotential_half_m2_s2.back(), 123.0, 1.0e-12);
  MPS_CHECK_NEAR(hydrostatic.residual_m2_s2[0], 0.0, 1.0e-10);
  MPS_CHECK_NEAR(hydrostatic.residual_m2_s2[1], 0.0, 1.0e-10);
}

MPS_TEST_CASE("isothermal and dry-adiabatic hydrostatic integrals are analytic") {
  const auto geometry = coordinate().geometry(100000.0, 10.0, 287.0, 1004.0, 100000.0);
  std::vector<mps::Real> isothermal_theta(2);
  for (std::size_t k = 0; k < isothermal_theta.size(); ++k) {
    isothermal_theta[k] = 280.0 / geometry.exner_full[k];
  }
  const auto isothermal =
      mps::integrate_hydrostatic_column(geometry, isothermal_theta, 1004.0, 10.0, 0.0);
  const std::vector<mps::Real> adiabatic_theta(2, 300.0);
  const auto adiabatic =
      mps::integrate_hydrostatic_column(geometry, adiabatic_theta, 1004.0, 10.0, 0.0);
  for (std::size_t k = 0; k < 2; ++k) {
    const double isothermal_exact =
        287.0 * 280.0 *
        std::log(geometry.pressure_half_pa[k + 1] / geometry.pressure_half_pa[k]);
    MPS_CHECK_NEAR(isothermal.geopotential_half_m2_s2[k] -
                       isothermal.geopotential_half_m2_s2[k + 1],
                   isothermal_exact, 2.0e-10 * isothermal_exact);
    const double adiabatic_exact =
        1004.0 * 300.0 * (geometry.exner_half[k + 1] - geometry.exner_half[k]);
    MPS_CHECK_NEAR(
        adiabatic.geopotential_half_m2_s2[k] - adiabatic.geopotential_half_m2_s2[k + 1],
        adiabatic_exact, 1.0e-10);
  }
}

MPS_TEST_CASE(
    "mass-flux recurrence closes and donor transport preserves constant scalars") {
  const std::vector<mps::Real> h{1.0, -0.25};
  const auto flux = mps::diagnose_vertical_mass_flux(h, {0.0, 0.5, 1.0}, 10.0);
  MPS_CHECK_NEAR(flux.interface_flux_kg_m2_s.front(), 0.0, 1.0e-15);
  MPS_CHECK_NEAR(flux.interface_flux_kg_m2_s.back(), 0.0, 1.0e-15);
  const std::vector<mps::Real> mass{100.0, 100.0};
  const std::vector<mps::Real> scalar_tendency{300.0, -75.0};
  const std::vector<mps::Real> scalar_mass{30000.0, 30000.0};
  const auto tendency = mps::vertical_scalar_rhs(
      scalar_mass, mass, scalar_tendency, flux,
      mps::VerticalTransportScheme::kDonorCell, mps::VerticalLimiterKind::kNone);
  MPS_CHECK_NEAR(tendency[0], 300.0 * flux.target_air_mass_tendency_kg_m2_s[0],
                 1.0e-12);
  MPS_CHECK_NEAR(tendency[1], 300.0 * flux.target_air_mass_tendency_kg_m2_s[1],
                 1.0e-12);
  MPS_CHECK_THROWS_AS(mps::diagnose_vertical_mass_flux(h, {0.0, 0.5, 0.9}, 10.0),
                      std::runtime_error);
}

MPS_TEST_CASE("donor and unlimited-linear transport reach their designed orders") {
  const double donor_coarse =
      spatial_transport_error(32, mps::VerticalTransportScheme::kDonorCell);
  const double donor_fine =
      spatial_transport_error(64, mps::VerticalTransportScheme::kDonorCell);
  const double linear_coarse =
      spatial_transport_error(32, mps::VerticalTransportScheme::kLinear);
  const double linear_fine =
      spatial_transport_error(64, mps::VerticalTransportScheme::kLinear);
  const double donor_order = std::log2(donor_coarse / donor_fine);
  const double linear_order = std::log2(linear_coarse / linear_fine);
  std::cout << "vertical transport L2 errors donor=" << donor_coarse << ','
            << donor_fine << " linear=" << linear_coarse << ',' << linear_fine
            << " orders=" << donor_order << ',' << linear_order << '\n';
  MPS_CHECK(donor_order >= 0.8);
  MPS_CHECK(linear_order >= 1.7);
}

MPS_TEST_CASE("minmod transport is conservative bounded and nonnegative") {
  const std::vector<mps::Real> air_mass{1.0, 2.0, 1.0, 3.0, 1.0, 2.0};
  const std::vector<mps::Real> scalar{0.0, 1.0, 0.25, 2.0, 0.5, 1.5};
  std::vector<mps::Real> scalar_mass(scalar.size());
  for (std::size_t k = 0; k < scalar.size(); ++k) {
    scalar_mass[k] = air_mass[k] * scalar[k];
  }
  mps::VerticalMassFlux flux;
  flux.interface_flux_kg_m2_s = {0.0, -0.2, -0.1, 0.15, 0.25, -0.1, 0.0};
  const std::vector<mps::Real> horizontal_tendency(scalar.size(), 0.0);
  const double dt = mps::vertical_stable_time_step(
      air_mass, flux.interface_flux_kg_m2_s, horizontal_tendency, 0.5, 10.0);
  const auto scalar_rhs = mps::vertical_scalar_rhs(
      scalar_mass, air_mass, horizontal_tendency, flux,
      mps::VerticalTransportScheme::kLinear, mps::VerticalLimiterKind::kMinmod);
  MPS_CHECK_NEAR(std::accumulate(scalar_rhs.begin(), scalar_rhs.end(), 0.0), 0.0,
                 1.0e-15);
  for (std::size_t k = 0; k < scalar.size(); ++k) {
    const double new_mass = air_mass[k] + dt * (flux.interface_flux_kg_m2_s[k] -
                                                flux.interface_flux_kg_m2_s[k + 1]);
    const double new_scalar = (scalar_mass[k] + dt * scalar_rhs[k]) / new_mass;
    MPS_CHECK(new_scalar >= 0.0);
    MPS_CHECK(new_scalar <= 2.0);
  }
}

MPS_TEST_CASE(
    "moving surface-pressure column preserves constant scalar and restart state") {
  const auto config = column_config();
  const auto coordinate_value = mps::make_vertical_coordinate(config);
  const auto initial =
      mps::make_vertical_column_initial_state(config, coordinate_value);
  const auto stopped = mps::run_vertical_column(config, initial, 20);
  MPS_CHECK(!stopped.reached_end_time);
  const auto restarted = mps::run_vertical_column(config, stopped.state);
  const auto direct = mps::run_vertical_column(config, initial);
  MPS_CHECK(restarted.reached_end_time);
  MPS_CHECK_EQ(restarted.state.surface_pressure_pa, direct.state.surface_pressure_pa);
  MPS_CHECK(restarted.state.potential_temperature_mass_k_kg_m2 ==
            direct.state.potential_temperature_mass_k_kg_m2);
  MPS_CHECK(restarted.state.tracer_mass_kg_m2 == direct.state.tracer_mass_kg_m2);
  MPS_CHECK_NEAR(direct.state.surface_pressure_pa, initial.surface_pressure_pa, 1.0e-9);
  const auto final_geometry = coordinate_value.geometry(
      direct.state.surface_pressure_pa, config.planet.gravity_m_s2,
      config.planet.gas_constant_j_kg_k, config.planet.heat_capacity_cp_j_kg_k,
      config.planet.reference_pressure_pa);
  for (std::size_t k = 0; k < direct.state.tracer_mass_kg_m2.size(); ++k) {
    MPS_CHECK_NEAR(direct.state.tracer_mass_kg_m2[k] / final_geometry.air_mass_kg_m2[k],
                   0.1, 256.0 * std::numeric_limits<double>::epsilon());
    MPS_CHECK_NEAR(direct.state.potential_temperature_mass_k_kg_m2[k] /
                       final_geometry.air_mass_kg_m2[k],
                   300.0, 256.0 * std::numeric_limits<double>::epsilon() * 300.0);
    MPS_CHECK_NEAR(direct.state.tracer_mass_kg_m2[k] / initial.tracer_mass_kg_m2[k],
                   direct.state.potential_temperature_mass_k_kg_m2[k] /
                       initial.potential_temperature_mass_k_kg_m2[k],
                   1.0e-10);
  }
}

MPS_TEST_CASE("driver honors start time and SSP-RK3 reaches temporal order") {
  auto config = column_config();
  config.run.start_time_s = 10.0;
  config.run.end_time_s = 30.0;
  const auto result = mps::run_vertical_column(config);
  MPS_CHECK_EQ(result.state.time_s, 30.0);
  MPS_CHECK_EQ(result.state.step, 80U);
  MPS_CHECK_NEAR(result.state.surface_pressure_pa, config.vertical.surface_pressure_pa,
                 1.0e-9);

  const double coarse = surface_pressure_quarter_period_error(2.0);
  const double medium = surface_pressure_quarter_period_error(1.0);
  const double fine = surface_pressure_quarter_period_error(0.5);
  const double coarse_order = std::log2(coarse / medium);
  const double fine_order = std::log2(medium / fine);
  std::cout << "SSP-RK3 surface-pressure errors=" << coarse << ',' << medium << ','
            << fine << " orders=" << coarse_order << ',' << fine_order << '\n';
  MPS_CHECK(coarse_order >= 2.7);
  MPS_CHECK(fine_order >= 2.7);
}

MPS_TEST_CASE("driver reduces a pressure-bound stage and reports physical CFL") {
  auto config = column_config();
  config.run.end_time_s = 60.0;
  config.run.time_step_s = 10.0;
  config.vertical.maximum_surface_pressure_pa = 100096.0;
  config.vertical.forcing_amplitude = 10.0;
  const auto result = mps::run_vertical_column(config);
  MPS_CHECK(result.reached_end_time);
  MPS_CHECK(result.maximum_cfl < 0.01);
  MPS_CHECK(result.maximum_cfl > 0.0);
  const auto diagnostics = mps::diagnostics::diagnose_vertical_column(
      config, mps::make_vertical_coordinate(config), result.state, result.mass_flux,
      result.maximum_cfl, result.budget);
  MPS_CHECK_NEAR(diagnostics.dry_mass_budget_residual_kg_m2, 0.0, 1.0e-8);
  MPS_CHECK_NEAR(diagnostics.potential_temperature_mass_budget_residual_k_kg_m2, 0.0,
                 1.0e-6);
  MPS_CHECK_NEAR(diagnostics.tracer_mass_budget_residual_kg_m2, 0.0, 1.0e-8);
  MPS_CHECK_NEAR(diagnostics.top_mass_flux_kg_m2_s, 0.0, 1.0e-15);
  MPS_CHECK_NEAR(diagnostics.surface_mass_flux_kg_m2_s, 0.0, 1.0e-15);
  MPS_CHECK_NEAR(diagnostics.continuity_residual_pa_s, 0.0, 1.0e-14);
  MPS_CHECK_EQ(diagnostics.non_finite_count, 0U);
}

MPS_TEST_CASE("static columns retain their prognostic values exactly") {
  auto config = column_config();
  config.vertical.test_case = mps::VerticalTestCase::kDryAdiabatic;
  config.vertical.forcing_amplitude = 0.0;
  const auto coordinate_value = mps::make_vertical_coordinate(config);
  const auto initial =
      mps::make_vertical_column_initial_state(config, coordinate_value);
  const auto result = mps::run_vertical_column(config, initial);
  MPS_CHECK_EQ(result.state.surface_pressure_pa, initial.surface_pressure_pa);
  MPS_CHECK_EQ(result.state.potential_temperature_mass_k_kg_m2[0],
               initial.potential_temperature_mass_k_kg_m2[0]);
  MPS_CHECK_EQ(result.state.tracer_mass_kg_m2[1], initial.tracer_mass_kg_m2[1]);
  const auto flat = mps::flatten_vertical_column_state(result.state);
  const auto restored = mps::unflatten_vertical_column_state(
      result.state.time_s, result.state.step, flat, coordinate_value.levels());
  MPS_CHECK_EQ(restored.tracer_mass_kg_m2[0], result.state.tracer_mass_kg_m2[0]);
}

MPS_TEST_CASE("vertical checkpoint round trips and rejects invalid restart time") {
  const auto config = column_config();
  const auto stopped = mps::run_vertical_column(config, std::nullopt, 20);
  const auto fingerprint = mps::config_fingerprint(config);
  std::stringstream stream;
  mps::write_checkpoint(
      stream, {.time_s = stopped.state.time_s,
               .step = stopped.state.step,
               .state = mps::flatten_vertical_column_state(stopped.state),
               .config_fingerprint = fingerprint,
               .layout_id = std::string(mps::kVerticalColumnCheckpointLayout)});
  const auto checkpoint =
      mps::read_checkpoint(stream, fingerprint, mps::kVerticalColumnCheckpointLayout,
                           1 + 2 * static_cast<std::size_t>(config.vertical.levels));
  const auto restored = mps::unflatten_vertical_column_state(
      checkpoint.time_s, checkpoint.step, checkpoint.state,
      static_cast<std::size_t>(config.vertical.levels));
  MPS_CHECK_EQ(restored.surface_pressure_pa, stopped.state.surface_pressure_pa);
  MPS_CHECK(restored.potential_temperature_mass_k_kg_m2 ==
            stopped.state.potential_temperature_mass_k_kg_m2);
  MPS_CHECK(restored.tracer_mass_kg_m2 == stopped.state.tracer_mass_kg_m2);

  auto invalid = restored;
  invalid.time_s = config.run.start_time_s - 1.0;
  MPS_CHECK_THROWS_AS(mps::run_vertical_column(config, invalid), std::invalid_argument);
}

MPS_TEST_CASE("uniform sigma coefficients are the only refinable coordinate family") {
  // ADR 0007: the ramp preserves the declared model top and satisfies the ADR 0005
  // endpoints exactly, not merely to within the rounding of the interior formula.
  const auto eight = mps::uniform_sigma_coefficients(1000.0, 8);
  MPS_CHECK_EQ(eight.a_half_pa.size(), 9U);
  MPS_CHECK_EQ(eight.a_half_pa.front(), 1000.0);
  MPS_CHECK_EQ(eight.a_half_pa.back(), 0.0);
  MPS_CHECK_EQ(eight.b_half.front(), 0.0);
  MPS_CHECK_EQ(eight.b_half.back(), 1.0);
  MPS_CHECK_EQ(eight.a_half_pa[1], 875.0);
  MPS_CHECK_EQ(eight.b_half[1], 0.125);
  // The shipped visualizer preset is exactly this ramp, so refining it reproduces the
  // preset when the requested level count equals the configured one.
  MPS_CHECK(mps::is_uniform_sigma(eight.a_half_pa, eight.b_half));
  const auto thirty = mps::uniform_sigma_coefficients(1000.0, 30);
  MPS_CHECK_EQ(thirty.a_half_pa.size(), 31U);
  MPS_CHECK(mps::is_uniform_sigma(thirty.a_half_pa, thirty.b_half));
  mps::HybridPressureCoefficients{thirty.a_half_pa, thirty.b_half}.validate(
      80000.0, 120000.0, 100.0);

  // A stretched coordinate is rejected as refinable even though it is a valid
  // coordinate.
  const std::vector<mps::Real> stretched_a{1000.0, 900.0, 500.0, 0.0};
  const std::vector<mps::Real> stretched_b{0.0, 0.1, 0.5, 1.0};
  MPS_CHECK(!mps::is_uniform_sigma(stretched_a, stretched_b));
  const std::vector<mps::Real> nudged_a{1000.0, 875.0, 750.0, 625.0, 500.0,
                                        375.0,  250.0, 125.0, 0.0};
  std::vector<mps::Real> nudged_b{0.0,   0.125, 0.25,  0.375, 0.5,
                                  0.625, 0.75,  0.875, 1.0};
  nudged_b[4] += 1.0e-6;
  MPS_CHECK(!mps::is_uniform_sigma(nudged_a, nudged_b));
  MPS_CHECK(!mps::is_uniform_sigma(std::vector<mps::Real>{1000.0},
                                   std::vector<mps::Real>{0.0}));
  MPS_CHECK_THROWS_AS(mps::uniform_sigma_coefficients(0.0, 8), std::invalid_argument);
  MPS_CHECK_THROWS_AS(mps::uniform_sigma_coefficients(1000.0, 0),
                      std::invalid_argument);
}

int main() { return mps::test::run_all(); }
