#include <cmath>
#include <numeric>
#include <vector>

#include "myplanetsim/dynamics/vertical_column_driver.hpp"
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
  MPS_CHECK_NEAR(restarted.state.surface_pressure_pa, direct.state.surface_pressure_pa,
                 1.0e-9);
  for (std::size_t k = 0; k < direct.state.tracer_mass_kg_m2.size(); ++k) {
    MPS_CHECK_NEAR(restarted.state.tracer_mass_kg_m2[k],
                   direct.state.tracer_mass_kg_m2[k], 1.0e-9);
    MPS_CHECK_NEAR(direct.state.tracer_mass_kg_m2[k] / initial.tracer_mass_kg_m2[k],
                   direct.state.potential_temperature_mass_k_kg_m2[k] /
                       initial.potential_temperature_mass_k_kg_m2[k],
                   1.0e-10);
  }
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

int main() { return mps::test::run_all(); }
