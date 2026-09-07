#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "myplanetsim/physics/moist_thermodynamics.hpp"
#include "myplanetsim/physics/saturation_adjustment.hpp"
#include "support/test.hpp"

namespace {

[[nodiscard]] mps::DiluteMoistThermodynamics thermodynamics() {
  return {.gas_constant_dry_air_j_kg_k = 287.0, .heat_capacity_cp_j_kg_k = 1004.0};
}

[[nodiscard]] double enthalpy(const std::vector<double>& mass,
                              const std::vector<double>& exner,
                              const std::vector<double>& theta_mass,
                              const std::vector<double>& vapor_mass) {
  const auto moist = thermodynamics();
  double result = 0.0;
  for (std::size_t level = 0; level < mass.size(); ++level)
    result += moist.heat_capacity_cp_j_kg_k * exner[level] * theta_mass[level] +
              moist.latent_heat_vaporization_j_kg * vapor_mass[level];
  return result;
}

}  // namespace

MPS_TEST_CASE("liquid saturation formula and analytic derivative agree") {
  const auto moist = thermodynamics();
  MPS_CHECK_NEAR(mps::saturation_vapor_pressure_pa(273.16, moist), 611.657, 1e-12);
  const double temperature = 285.0;
  const double pressure = 85000.0;
  const double step = 1e-3;
  const double finite_difference =
      (mps::saturation_mixing_ratio(temperature + step, pressure, moist) -
       mps::saturation_mixing_ratio(temperature - step, pressure, moist)) /
      (2.0 * step);
  const double analytic = mps::saturation_mixing_ratio_temperature_derivative_k_1(
      temperature, pressure, moist);
  MPS_CHECK_NEAR(analytic, finite_difference, 1e-8 * analytic);
  MPS_CHECK_NEAR(mps::relative_humidity(
                     0.6 * mps::saturation_mixing_ratio(temperature, pressure, moist),
                     temperature, pressure, moist),
                 0.6, 1e-14);
}

MPS_TEST_CASE("dry subsaturated and exactly saturated layers are unchanged") {
  const auto moist = thermodynamics();
  std::vector<double> mass{10.0, 20.0, 30.0};
  std::vector<double> pressure{30000.0, 60000.0, 90000.0};
  std::vector<double> exner{0.7, 0.85, 0.97};
  std::vector<double> theta_mass{mass[0] * 250.0 / exner[0], mass[1] * 270.0 / exner[1],
                                 mass[2] * 290.0 / exner[2]};
  std::vector<double> vapor_mass{
      0.0, 0.5 * mass[1] * mps::saturation_mixing_ratio(270.0, pressure[1], moist),
      mass[2] * mps::saturation_mixing_ratio(290.0, pressure[2], moist)};
  const auto theta_before = theta_mass;
  const auto vapor_before = vapor_mass;
  const auto result = mps::adjust_saturation_column(mass, pressure, exner, theta_mass,
                                                    vapor_mass, moist);
  MPS_CHECK_EQ(result.adjusted_layer_count, 0U);
  MPS_CHECK_EQ(result.condensed_water_kg_m2, 0.0);
  MPS_CHECK(theta_mass == theta_before);
  MPS_CHECK(vapor_mass == vapor_before);
}

MPS_TEST_CASE("supersaturation warms to saturation and conserves enthalpy") {
  const auto moist = thermodynamics();
  std::vector<double> mass{8.0, 23.0, 51.0};
  std::vector<double> pressure{35000.0, 65000.0, 95000.0};
  std::vector<double> exner{0.72, 0.88, 0.98};
  std::vector<double> temperature{245.0, 275.0, 292.0};
  std::vector<double> theta_mass(3);
  std::vector<double> vapor_mass(3);
  for (std::size_t level = 0; level < 3; ++level) {
    theta_mass[level] = mass[level] * temperature[level] / exner[level];
    vapor_mass[level] = mass[level] * (mps::saturation_mixing_ratio(
                                           temperature[level], pressure[level], moist) +
                                       0.001 * static_cast<double>(level + 1));
  }
  const double before = enthalpy(mass, exner, theta_mass, vapor_mass);
  const auto result = mps::adjust_saturation_column(mass, pressure, exner, theta_mass,
                                                    vapor_mass, moist);
  const double after = enthalpy(mass, exner, theta_mass, vapor_mass);
  MPS_CHECK_EQ(result.adjusted_layer_count, 3U);
  MPS_CHECK(result.condensed_water_kg_m2 > 0.0);
  MPS_CHECK_NEAR(after, before, 1e-12 * before);
  MPS_CHECK_NEAR(result.enthalpy_change_j_m2, after - before, 1e-13 * before);
  for (std::size_t level = 0; level < 3; ++level) {
    const double final_temperature = theta_mass[level] / mass[level] * exner[level];
    const double final_vapor = vapor_mass[level] / mass[level];
    const double saturation =
        mps::saturation_mixing_ratio(final_temperature, pressure[level], moist);
    MPS_CHECK(final_vapor >= 0.0);
    MPS_CHECK(final_vapor - saturation <= 1e-12 + 1e-10 * saturation);
  }

  const auto theta_once = theta_mass;
  const auto vapor_once = vapor_mass;
  const auto second = mps::adjust_saturation_column(mass, pressure, exner, theta_mass,
                                                    vapor_mass, moist);
  MPS_CHECK_EQ(second.adjusted_layer_count, 0U);
  MPS_CHECK(theta_mass == theta_once);
  MPS_CHECK(vapor_mass == vapor_once);
}

MPS_TEST_CASE("invalid and non-dilute inputs are rejected without clipping") {
  const auto moist = thermodynamics();
  MPS_CHECK_THROWS_AS(mps::validate_dilute_moist_state(300.0, 100000.0, 0.11, moist),
                      std::domain_error);
  MPS_CHECK_THROWS_AS(mps::validate_dilute_moist_state(0.0, 100000.0, 0.01, moist),
                      std::invalid_argument);
}

int main() { return mps::test::run_all(); }
