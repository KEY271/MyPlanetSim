#include <cmath>
#include <sstream>

#include "myplanetsim/diagnostics/planetary_scales.hpp"
#include "support/test.hpp"

MPS_TEST_CASE("planetary scale definitions use explicit input scales") {
  const auto planet = mps::PlanetParameters::earth_like();
  const mps::PlanetaryScaleInputs inputs{
      .characteristic_wind_m_s = 20.0,
      .characteristic_length_m = 1.0e6,
      .buoyancy_frequency_s_1 = 0.01,
      .radiative_time_s = 40.0 * 86400.0,
      .reference_temperature_k = 288.0,
      .temperature_contrast_k = 60.0,
  };
  const auto scales = mps::compute_planetary_scales(planet, inputs);
  const auto f_reference = std::sqrt(2.0) * planet.rotation_rate_rad_s;
  const auto height = planet.gas_constant_j_kg_k * 288.0 / planet.gravity_m_s2;
  MPS_CHECK(scales.rossby_number.defined);
  MPS_CHECK_NEAR(scales.rossby_number.value, 20.0 / (f_reference * 1.0e6), 1.0e-14);
  MPS_CHECK_NEAR(scales.deformation_radius_m.value, 0.01 * height / f_reference,
                 1.0e-9);
  MPS_CHECK(scales.burger_number.defined);
  MPS_CHECK(scales.thermal_rossby_number.defined);
  MPS_CHECK(scales.radiative_over_rotation.defined);
}

MPS_TEST_CASE("zero rotation uses explicit undefined values") {
  auto planet = mps::PlanetParameters::earth_like();
  planet.rotation_rate_rad_s = 0.0;
  const mps::PlanetaryScaleInputs inputs{20.0, 1.0e6, 0.01, 86400.0, 288.0, 60.0};
  const auto scales = mps::compute_planetary_scales(planet, inputs);
  MPS_CHECK(!scales.rossby_number.defined);
  MPS_CHECK(!scales.deformation_radius_m.defined);
  MPS_CHECK_EQ(scales.rossby_number.value, 0.0);
  MPS_CHECK(std::isfinite(scales.rossby_number.value));
}

MPS_TEST_CASE("zero emissivity has undefined surface radiative time") {
  const auto undefined = mps::surface_radiative_time_scale(2.0e6, 0.0, 288.0);
  MPS_CHECK(!undefined.defined);
  const auto land = mps::surface_radiative_time_scale(2.0e6, 1.0, 288.0);
  const auto ocean = mps::surface_radiative_time_scale(4.0e7, 1.0, 288.0);
  MPS_CHECK(land.defined);
  MPS_CHECK_NEAR(ocean.value / land.value, 20.0, 1.0e-14);
}

MPS_TEST_CASE("undefined metadata omits a fabricated numeric value") {
  std::ostringstream output;
  mps::write_defined_quantity(output, "planetary.rossby_number", {});
  MPS_CHECK_EQ(output.str(), "planetary.rossby_number.defined = false\n");
}

int main() { return mps::test::run_all(); }
