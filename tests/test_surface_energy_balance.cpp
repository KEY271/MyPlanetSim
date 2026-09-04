#include <cmath>
#include <vector>

#include "myplanetsim/dynamics/dry_hydrostatic_driver.hpp"
#include "myplanetsim/physics/surface_energy_balance.hpp"
#include "support/test.hpp"

namespace {

struct Fixture {
  mps::CubedSphereGrid grid{1, 2.0};
  mps::DryHydrostaticDerived atmosphere;
  std::vector<mps::Real> surface_pressure;
  mps::OrbitState orbit{{1.0, 0.0, 0.0}, 1000.0, 1.0, 0.0};

  Fixture() {
    atmosphere.cells = grid.cell_count();
    atmosphere.levels = 1;
    atmosphere.pressure_pa.assign(grid.cell_count(), 100000.0);
    atmosphere.air_mass_kg_m2.assign(grid.cell_count(), 1000.0);
    atmosphere.velocity_m_s.assign(grid.cell_count(), {});
    atmosphere.temperature_k.assign(grid.cell_count(), 280.0);
    surface_pressure.assign(grid.cell_count(), 100000.0);
  }
};

[[nodiscard]] mps::SurfaceParameters parameters() {
  return {.geography = mps::SurfaceGeography::kUniform,
          .uniform_land_fraction = 1.0,
          .land_heat_capacity_j_m2_k = 2.0e6,
          .ocean_heat_capacity_j_m2_k = 4.0e7,
          .initial_temperature_k = 300.0,
          .albedo = 0.3,
          .emissivity = 1.0,
          .air_exchange_coefficient_w_m2_k = 10.0,
          .internal_heat_flux_w_m2 = 0.0};
}

[[nodiscard]] mps::SurfaceBoundary boundary(const Fixture& fixture,
                                            const mps::Real land_fraction) {
  return mps::SurfaceBoundary(
      std::vector<mps::Real>(fixture.grid.cell_count(), 0.0),
      std::vector<mps::Real>(fixture.grid.cell_count(), land_fraction), "test",
      "test:fingerprint");
}

}  // namespace

MPS_TEST_CASE("surface energy signs and budget close from one evaluated flux") {
  const Fixture fixture;
  const auto surface = boundary(fixture, 1.0);
  const std::vector<mps::Real> temperature(fixture.grid.cell_count(), 300.0);
  const auto tendency = mps::surface_energy_tendency(
      fixture.grid, surface, temperature, fixture.atmosphere, fixture.surface_pressure,
      mps::PlanetParameters::earth_like(), parameters(), fixture.orbit);
  MPS_CHECK(tendency.diagnostics.absorbed_stellar_power_w > 0.0);
  MPS_CHECK(tendency.diagnostics.outgoing_longwave_power_w > 0.0);
  MPS_CHECK(tendency.diagnostics.sensible_to_atmosphere_power_w > 0.0);
  MPS_CHECK_NEAR(tendency.diagnostics.surface_budget_residual_w, 0.0, 1.0e-8);
  for (const auto source : tendency.potential_temperature_mass_k_kg_m2_s)
    MPS_CHECK(source > 0.0);
}

MPS_TEST_CASE("ocean responds less than land to the same forcing") {
  const Fixture fixture;
  const std::vector<mps::Real> temperature(fixture.grid.cell_count(), 300.0);
  const auto land = mps::surface_energy_tendency(
      fixture.grid, boundary(fixture, 1.0), temperature, fixture.atmosphere,
      fixture.surface_pressure, mps::PlanetParameters::earth_like(), parameters(),
      fixture.orbit);
  const auto ocean = mps::surface_energy_tendency(
      fixture.grid, boundary(fixture, 0.0), temperature, fixture.atmosphere,
      fixture.surface_pressure, mps::PlanetParameters::earth_like(), parameters(),
      fixture.orbit);
  MPS_CHECK_NEAR(
      land.surface_temperature_k_s.front() / ocean.surface_temperature_k_s.front(),
      20.0, 1.0e-12);
}

MPS_TEST_CASE("no-insolation radiative equilibrium is exact") {
  const Fixture fixture;
  auto configured = parameters();
  configured.albedo = 1.0;
  configured.air_exchange_coefficient_w_m2_k = 0.0;
  configured.internal_heat_flux_w_m2 = mps::kStefanBoltzmannWm2K4 * std::pow(300.0, 4);
  const std::vector<mps::Real> temperature(fixture.grid.cell_count(), 300.0);
  const auto tendency = mps::surface_energy_tendency(
      fixture.grid, boundary(fixture, 1.0), temperature, fixture.atmosphere,
      fixture.surface_pressure, mps::PlanetParameters::earth_like(), configured,
      fixture.orbit);
  for (const auto rate : tendency.surface_temperature_k_s) MPS_CHECK_EQ(rate, 0.0);
}

MPS_TEST_CASE("driver advances the coupled surface at SSP-RK3 stages") {
  mps::ExperimentConfig config{
      .kind = mps::ExperimentKind::kDryHydrostatic,
      .planet = {6371220, 7.29212e-5, 9.80616, 287, 1004, 100000},
      .run = {0, 2, 1, 0},
      .grid = {2},
      .vertical = {.levels = 2,
                   .a_half_pa = {1000, 500, 0},
                   .b_half = {0, .5, 1},
                   .surface_pressure_pa = 100000,
                   .minimum_surface_pressure_pa = 80000,
                   .maximum_surface_pressure_pa = 120000,
                   .minimum_pressure_thickness_pa = 100,
                   .initial_temperature_k = 280,
                   .initial_potential_temperature_k = 300,
                   .temperature_floor_k = 100,
                   .transport_scheme = mps::VerticalTransportScheme::kLinear,
                   .limiter = mps::VerticalLimiterKind::kMinmod,
                   .cfl = .5},
      .dry_hydrostatic = {.test_case = mps::DryHydrostaticTestCase::kIsothermalRest},
      .physics = {.kind = mps::PhysicsKind::kSurfaceEnergyBalance},
      .diagnostics = {1},
      .output_directory = "x"};
  config.orbit = mps::OrbitParameters{86400, 0, 0, 0, 0, 0, 1361};
  config.surface = mps::SurfaceParameters{.geography = mps::SurfaceGeography::kUniform,
                                          .uniform_land_fraction = 0.5,
                                          .land_heat_capacity_j_m2_k = 2e6,
                                          .ocean_heat_capacity_j_m2_k = 4e7,
                                          .initial_temperature_k = 288,
                                          .albedo = 0.3,
                                          .emissivity = 1,
                                          .air_exchange_coefficient_w_m2_k = 10,
                                          .internal_heat_flux_w_m2 = 0};
  config.validate();
  const mps::DryHydrostaticDriver driver(config);
  auto state = driver.initial_state();
  const auto initial_surface = state.surface_temperature_k;
  const auto initial_atmosphere = state.potential_temperature_mass_k_kg_m2;
  driver.advance(state, 1.0);
  MPS_CHECK_EQ(state.surface_temperature_k.size(), driver.grid().cell_count());
  MPS_CHECK(state.surface_temperature_k != initial_surface);
  MPS_CHECK(state.potential_temperature_mass_k_kg_m2 != initial_atmosphere);
  for (const auto temperature : state.surface_temperature_k)
    MPS_CHECK(temperature > 0.0);
}

// ADR 0011: the reservoir is integrated explicitly, so its step is bounded before the
// fact. Previously only positivity was checked, and a small heat capacity could
// oscillate inside the positive range without failing any gate.
MPS_TEST_CASE("the surface reservoir reports its own stability limit") {
  const Fixture fixture;
  auto small = parameters();
  small.land_heat_capacity_j_m2_k = 1.0e4;
  const std::vector<mps::Real> temperature(fixture.grid.cell_count(), 300.0);
  const auto tendency = mps::surface_energy_tendency(
      fixture.grid, boundary(fixture, 1.0), temperature, fixture.atmosphere,
      fixture.surface_pressure, {2.0, 1.0e-4, 9.8, 287.0, 1004.0, 100000.0}, small,
      fixture.orbit);

  const mps::Real relaxation =
      4.0 * small.emissivity * mps::kStefanBoltzmannWm2K4 * std::pow(300.0, 3);
  MPS_CHECK_NEAR(tendency.stable_time_step_s,
                 small.cfl * small.land_heat_capacity_j_m2_k / relaxation, 1.0e-9);

  // A twenty-times larger reservoir relaxes twenty times more slowly.
  auto large = small;
  large.land_heat_capacity_j_m2_k = 2.0e5;
  const auto slower = mps::surface_energy_tendency(
      fixture.grid, boundary(fixture, 1.0), temperature, fixture.atmosphere,
      fixture.surface_pressure, {2.0, 1.0e-4, 9.8, 287.0, 1004.0, 100000.0}, large,
      fixture.orbit);
  MPS_CHECK_NEAR(slower.stable_time_step_s, 20.0 * tendency.stable_time_step_s, 1.0e-6);
}

int main() { return mps::test::run_all(); }
