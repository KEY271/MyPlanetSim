#include <cmath>

#include "myplanetsim/dynamics/dry_hydrostatic_driver.hpp"
#include "myplanetsim/physics/planetary_newtonian.hpp"
#include "support/test.hpp"

namespace {

[[nodiscard]] mps::ExperimentConfig driver_config(const mps::PhysicsKind physics) {
  mps::ExperimentConfig config{
      .kind = mps::ExperimentKind::kDryHydrostatic,
      .planet = {6371220, 7.29212e-5, 9.80616, 287, 1004, 100000},
      .run = {0, 10, 1, 0},
      .grid = {2},
      .vertical = {.levels = 2,
                   .a_half_pa = {1000, 500, 0},
                   .b_half = {0, .5, 1},
                   .surface_pressure_pa = 100000,
                   .minimum_surface_pressure_pa = 80000,
                   .maximum_surface_pressure_pa = 120000,
                   .minimum_pressure_thickness_pa = 100,
                   .initial_temperature_k = 264,
                   .initial_potential_temperature_k = 300,
                   .temperature_floor_k = 100,
                   .transport_scheme = mps::VerticalTransportScheme::kLinear,
                   .limiter = mps::VerticalLimiterKind::kMinmod,
                   .cfl = .5},
      .dry_hydrostatic = {.test_case = mps::DryHydrostaticTestCase::kHeldSuarez},
      .physics = {.kind = physics},
      .diagnostics = {1},
      .output_directory = "x"};
  if (physics == mps::PhysicsKind::kPlanetaryNewtonian)
    config.surface = mps::SurfaceParameters{
        .geography = mps::SurfaceGeography::kUniform, .uniform_land_fraction = 0.0};
  return config;
}

}  // namespace

MPS_TEST_CASE("axisymmetric rates exactly preserve Held-Suarez") {
  const auto planet = mps::PlanetParameters::earth_like();
  for (const mps::Real latitude : {-1.0, -0.3, 0.0, 0.7, 1.2}) {
    for (const mps::Real pressure : {1000.0, 70000.0, 100000.0}) {
      const auto held = mps::held_suarez_rates(latitude, pressure, 100000.0, planet);
      const auto planetary = mps::planetary_newtonian_rates(
          latitude, -std::pow(std::sin(latitude), 2), pressure, 100000.0, planet);
      MPS_CHECK_EQ(planetary.equilibrium_temperature_k, held.equilibrium_temperature_k);
      MPS_CHECK_EQ(planetary.temperature_relaxation_rate_s_1,
                   held.temperature_relaxation_rate_s_1);
      MPS_CHECK_EQ(planetary.rayleigh_drag_rate_s_1, held.rayleigh_drag_rate_s_1);
    }
  }
}

MPS_TEST_CASE("substellar forcing retains the signed nightside dipole") {
  const auto planet = mps::PlanetParameters::earth_like();
  const auto day = mps::planetary_newtonian_rates(0.0, 1.0, 100000.0, 100000.0, planet);
  const auto night =
      mps::planetary_newtonian_rates(0.0, -1.0, 100000.0, 100000.0, planet);
  const auto clipped_night =
      mps::planetary_newtonian_rates(0.0, 0.0, 100000.0, 100000.0, planet);
  MPS_CHECK_EQ(day.equilibrium_temperature_k, 375.0);
  MPS_CHECK_EQ(night.equilibrium_temperature_k, 255.0);
  MPS_CHECK(night.equilibrium_temperature_k < clipped_night.equilibrium_temperature_k);
}

MPS_TEST_CASE("axisymmetric driver path is the Held-Suarez source path") {
  const mps::DryHydrostaticDriver held_driver(
      driver_config(mps::PhysicsKind::kHeldSuarez));
  const mps::DryHydrostaticDriver planetary_driver(
      driver_config(mps::PhysicsKind::kPlanetaryNewtonian));
  const auto held = held_driver.rhs(held_driver.initial_state());
  const auto planetary = planetary_driver.rhs(planetary_driver.initial_state());
  for (std::size_t index = 0; index < held.tendency.momentum.size(); ++index)
    MPS_CHECK_EQ(
        mps::norm(held.tendency.momentum[index] - planetary.tendency.momentum[index]),
        0.0);
  MPS_CHECK(held.tendency.potential_temperature_mass ==
            planetary.tendency.potential_temperature_mass);
  MPS_CHECK_EQ(held.physics_diagnostics.thermal_energy_rate_w,
               planetary.physics_diagnostics.thermal_energy_rate_w);
}

int main() { return mps::test::run_all(); }
