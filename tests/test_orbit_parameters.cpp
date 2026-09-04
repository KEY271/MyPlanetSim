#include <cmath>
#include <numbers>

#include "myplanetsim/core/orbit_parameters.hpp"
#include "myplanetsim/grid/cubed_sphere_grid.hpp"
#include "support/test.hpp"

namespace {

[[nodiscard]] mps::OrbitParameters circular_orbit() {
  return {.period_s = 100.0,
          .eccentricity = 0.0,
          .obliquity_rad = 0.0,
          .longitude_of_periapsis_rad = 0.0,
          .initial_mean_anomaly_rad = 0.0,
          .initial_substellar_longitude_rad = 0.4,
          .stellar_flux_at_semimajor_axis_w_m2 = 1000.0};
}

}  // namespace

MPS_TEST_CASE("initial substellar longitude anchors the rotating frame") {
  const auto state = mps::evaluate_orbit(circular_orbit(), 0.01, 0.0);
  MPS_CHECK_NEAR(std::atan2(state.star_direction_planet_fixed.y,
                            state.star_direction_planet_fixed.x),
                 0.4, 1.0e-14);
  MPS_CHECK_NEAR(state.star_direction_planet_fixed.z, 0.0, 1.0e-15);
  MPS_CHECK_NEAR(state.stellar_flux_w_m2, 1000.0, 1.0e-12);
}

MPS_TEST_CASE("substellar longitude moves at orbital minus rotation rate") {
  const auto orbit = circular_orbit();
  constexpr mps::Real rotation = 0.01;
  constexpr mps::Real time = 4.0;
  const auto state = mps::evaluate_orbit(orbit, rotation, time);
  const mps::Real expected = orbit.initial_substellar_longitude_rad +
                             (orbit.mean_motion_rad_s() - rotation) * time;
  MPS_CHECK_NEAR(std::atan2(state.star_direction_planet_fixed.y,
                            state.star_direction_planet_fixed.x),
                 expected, 2.0e-14);
}

MPS_TEST_CASE("synchronous circular orbit has a fixed star direction") {
  const auto orbit = circular_orbit();
  const auto rotation = orbit.mean_motion_rad_s();
  const auto initial = mps::evaluate_orbit(orbit, rotation, 0.0);
  const auto later = mps::evaluate_orbit(orbit, rotation, 37.0);
  MPS_CHECK_NEAR(mps::norm(later.star_direction_planet_fixed -
                           initial.star_direction_planet_fixed),
                 0.0, 2.0e-14);
  const auto day = mps::stellar_day(orbit, rotation);
  MPS_CHECK(day.synchronous);
  MPS_CHECK(!day.period_s.has_value());
}

MPS_TEST_CASE("eccentric flux ratio matches periapsis and apoapsis") {
  auto orbit = circular_orbit();
  orbit.eccentricity = 0.2;
  orbit.initial_mean_anomaly_rad = 0.0;
  const auto periapsis = mps::evaluate_orbit(orbit, 0.0, 0.0);
  const auto apoapsis = mps::evaluate_orbit(orbit, 0.0, 0.5 * orbit.period_s);
  MPS_CHECK_NEAR(periapsis.stellar_flux_w_m2 / apoapsis.stellar_flux_w_m2,
                 std::pow((1.0 + orbit.eccentricity) / (1.0 - orbit.eccentricity), 2),
                 2.0e-14);
}

MPS_TEST_CASE("solar zenith cosine clips the night side exactly") {
  const auto state = mps::evaluate_orbit(circular_orbit(), 0.0, 0.0);
  MPS_CHECK_NEAR(mps::cosine_solar_zenith(state.star_direction_planet_fixed, state),
                 1.0, 1.0e-15);
  MPS_CHECK_EQ(mps::cosine_solar_zenith(-state.star_direction_planet_fixed, state),
               0.0);
}

MPS_TEST_CASE("cubed sphere integrates instantaneous incoming power") {
  const auto orbit = circular_orbit();
  const auto state = mps::evaluate_orbit(orbit, 0.0, 0.0);
  constexpr mps::Real radius = 2.0;
  const mps::CubedSphereGrid grid(32, radius);
  mps::Real integrated = 0.0;
  for (const auto& cell : grid.cells()) {
    integrated += cell.area_m2 * mps::cosine_solar_zenith(cell.center, state) *
                  state.stellar_flux_w_m2;
  }
  const auto exact =
      std::numbers::pi_v<mps::Real> * radius * radius * state.stellar_flux_w_m2;
  MPS_CHECK_NEAR(integrated / exact, 1.0, 4.0e-4);
}

MPS_TEST_CASE("orbit rejects invalid eccentricity") {
  auto orbit = circular_orbit();
  orbit.eccentricity = 1.0;
  MPS_CHECK_THROWS_AS(orbit.validate(), std::invalid_argument);
}

int main() { return mps::test::run_all(); }
