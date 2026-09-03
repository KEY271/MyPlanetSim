#include <stdexcept>

#include "myplanetsim/core/planet_parameters.hpp"
#include "support/test.hpp"

MPS_TEST_CASE("Earth-like parameters are valid") {
  const auto planet = mps::PlanetParameters::earth_like();
  planet.validate();

  MPS_CHECK_NEAR(planet.radius_m, 6.37122e6, 1.0e-6);
  MPS_CHECK_NEAR(planet.heat_capacity_cv_j_kg_k(), 717.5, 1.0e-12);
  MPS_CHECK_NEAR(planet.kappa(), 287.0 / 1004.5, 1.0e-15);
  MPS_CHECK(planet.rotation_period_s().has_value());
  MPS_CHECK_NEAR(*planet.rotation_period_s(), 86164.1, 0.1);
}

MPS_TEST_CASE("planet validation rejects nonphysical thermodynamics") {
  auto planet = mps::PlanetParameters::earth_like();
  planet.radius_m = 0.0;
  MPS_CHECK_THROWS_AS(planet.validate(), std::invalid_argument);

  planet = mps::PlanetParameters::earth_like();
  planet.heat_capacity_cp_j_kg_k = planet.gas_constant_j_kg_k;
  MPS_CHECK_THROWS_AS(planet.validate(), std::invalid_argument);
}

MPS_TEST_CASE("zero and reverse rotation are supported") {
  auto planet = mps::PlanetParameters::earth_like();
  planet.rotation_rate_rad_s = 0.0;
  planet.validate();
  MPS_CHECK(!planet.rotation_period_s().has_value());

  planet.rotation_rate_rad_s = -7.292115e-5;
  planet.validate();
  MPS_CHECK(planet.rotation_period_s().has_value());
  MPS_CHECK_NEAR(*planet.rotation_period_s(), 86164.1, 0.1);
}

int main() { return mps::test::run_all(); }
