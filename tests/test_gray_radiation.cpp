#include <cmath>
#include <limits>
#include <vector>

#include "myplanetsim/physics/gray_radiation.hpp"
#include "support/test.hpp"

namespace {

[[nodiscard]] mps::RadiationParameters parameters() {
  return mps::RadiationParameters{
      .shortwave_absorption_m2_kg = 2.0e-5,
      .longwave_absorption_ref_m2_kg = 3.0e-4,
      .reference_pressure_pa = 100000.0,
      .longwave_pressure_exponent = 1.0,
      .longwave_diffusivity_factor = 1.66,
      .shortwave_diffuse_factor = 1.66,
      .cfl = 0.5,
  };
}

}  // namespace

MPS_TEST_CASE("gray optical depth follows mass path and telescopes") {
  constexpr mps::Real gravity = 10.0;
  const std::vector<mps::Real> pressure{1000.0, 10000.0, 40000.0, 101000.0};
  const auto opacity =
      mps::gray_radiation_optical_depth(pressure, gravity, parameters());
  MPS_CHECK_EQ(opacity.shortwave.size(), 3U);
  MPS_CHECK_EQ(opacity.longwave.size(), 3U);

  mps::Real shortwave_total = 0.0;
  mps::Real longwave_total = 0.0;
  for (std::size_t level = 0; level < opacity.shortwave.size(); ++level) {
    shortwave_total += opacity.shortwave[level];
    longwave_total += opacity.longwave[level];
  }
  const auto p = parameters();
  const mps::Real expected =
      p.shortwave_absorption_m2_kg * (pressure.back() - pressure.front()) / gravity;
  MPS_CHECK_NEAR(shortwave_total, expected, 1e-15);
  MPS_CHECK_NEAR(
      longwave_total,
      p.longwave_absorption_ref_m2_kg * (pressure.back() - pressure.front()) / gravity,
      1e-15);
}

MPS_TEST_CASE("gray longwave pressure exponent uses absolute pressure") {
  auto p = parameters();
  p.longwave_pressure_exponent = 2.0;
  constexpr mps::Real gravity = 9.81;
  const std::vector<mps::Real> low_surface{1000.0, 50000.0, 80000.0};
  const std::vector<mps::Real> high_surface{1000.0, 50000.0, 100000.0};
  const auto low = mps::gray_radiation_optical_depth(low_surface, gravity, p);
  const auto high = mps::gray_radiation_optical_depth(high_surface, gravity, p);
  const mps::Real low_total = low.longwave[0] + low.longwave[1];
  const mps::Real high_total = high.longwave[0] + high.longwave[1];
  const auto analytic = [&](const mps::Real bottom) {
    return p.longwave_absorption_ref_m2_kg * p.reference_pressure_pa / (2.0 * gravity) *
           (std::pow(bottom / p.reference_pressure_pa, 2.0) -
            std::pow(low_surface.front() / p.reference_pressure_pa, 2.0));
  };
  MPS_CHECK_NEAR(low_total, analytic(low_surface.back()), 1e-15);
  MPS_CHECK_NEAR(high_total, analytic(high_surface.back()), 1e-15);
  MPS_CHECK(high_total > low_total);
}

MPS_TEST_CASE("gray optical depth preserves a thin pressure layer") {
  auto p = parameters();
  p.longwave_pressure_exponent = 2.0;
  const std::vector<mps::Real> pressure{99999.999, 100000.0};
  const auto opacity = mps::gray_radiation_optical_depth(pressure, 10.0, p);
  const mps::Real expected =
      p.longwave_absorption_ref_m2_kg * (pressure.back() - pressure.front()) / 10.0;
  MPS_CHECK_NEAR(opacity.longwave.front(), expected, 1e-7 * expected);
  MPS_CHECK(opacity.longwave.front() > 0.0);
}

MPS_TEST_CASE("gray optical depth rejects invalid columns and parameters") {
  auto p = parameters();
  MPS_CHECK_THROWS_AS(
      mps::gray_radiation_optical_depth(std::vector<mps::Real>{1000.0}, 10.0, p),
      std::invalid_argument);
  MPS_CHECK_THROWS_AS(mps::gray_radiation_optical_depth(
                          std::vector<mps::Real>{1000.0, 1000.0}, 10.0, p),
                      std::invalid_argument);
  p.longwave_pressure_exponent = 0.5;
  MPS_CHECK_THROWS_AS(
      mps::gray_radiation_optical_depth(std::vector<mps::Real>{0.0, 100000.0}, 10.0, p),
      std::invalid_argument);
  p = parameters();
  p.shortwave_absorption_m2_kg = std::numeric_limits<mps::Real>::quiet_NaN();
  MPS_CHECK_THROWS_AS(
      mps::gray_radiation_optical_depth(std::vector<mps::Real>{0.0, 100000.0}, 10.0, p),
      std::invalid_argument);
}

int main() { return mps::test::run_all(); }
