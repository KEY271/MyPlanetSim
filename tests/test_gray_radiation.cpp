#include <cmath>
#include <limits>
#include <vector>

#include "myplanetsim/physics/gray_radiation.hpp"
#include "myplanetsim/physics/surface_energy_balance.hpp"
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

[[nodiscard]] mps::GrayRadiationColumnInput column_input(
    const std::vector<mps::Real>& pressure, const std::vector<mps::Real>& temperature,
    const mps::RadiationParameters& p) {
  return mps::GrayRadiationColumnInput{
      .pressure_half_pa = pressure,
      .temperature_k = temperature,
      .surface_temperature_k = 300.0,
      .gravity_m_s2 = 10.0,
      .stellar_flux_w_m2 = 1000.0,
      .cosine_solar_zenith = 0.5,
      .surface_albedo = 0.25,
      .surface_emissivity = 1.0,
      .parameters = p,
  };
}

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

MPS_TEST_CASE("transparent gray column preserves directional fluxes") {
  auto p = parameters();
  p.shortwave_absorption_m2_kg = 0.0;
  p.longwave_absorption_ref_m2_kg = 0.0;
  const std::vector<mps::Real> pressure{1000.0, 30000.0, 100000.0};
  const std::vector<mps::Real> temperature{220.0, 280.0};
  const auto result =
      mps::gray_radiation_column(column_input(pressure, temperature, p));
  const mps::Real surface_emission = mps::kStefanBoltzmannWm2K4 * std::pow(300.0, 4.0);
  for (std::size_t interface = 0; interface < pressure.size(); ++interface) {
    MPS_CHECK_NEAR(result.shortwave_down_w_m2[interface], 500.0, 0.0);
    MPS_CHECK_NEAR(result.shortwave_up_w_m2[interface], 125.0, 0.0);
    MPS_CHECK_NEAR(result.longwave_down_w_m2[interface], 0.0, 0.0);
    MPS_CHECK_NEAR(result.longwave_up_w_m2[interface], surface_emission, 1e-12);
  }
  for (const mps::Real convergence : result.radiative_convergence_w_m2)
    MPS_CHECK_NEAR(convergence, 0.0, 1e-12);
}

MPS_TEST_CASE("isothermal longwave sweep matches the analytic attenuation") {
  auto p = parameters();
  p.shortwave_absorption_m2_kg = 0.0;
  p.longwave_absorption_ref_m2_kg = 2e-4;
  p.longwave_pressure_exponent = 1.0;
  const std::vector<mps::Real> pressure{0.0, 25000.0, 60000.0, 100000.0};
  const std::vector<mps::Real> temperature{250.0, 250.0, 250.0};
  auto input = column_input(pressure, temperature, p);
  input.stellar_flux_w_m2 = 0.0;
  input.surface_temperature_k = 310.0;
  const auto result = mps::gray_radiation_column(input);
  const mps::Real source = mps::kStefanBoltzmannWm2K4 * std::pow(250.0, 4.0);
  const mps::Real boundary = mps::kStefanBoltzmannWm2K4 * std::pow(310.0, 4.0);

  mps::Real path_from_top = 0.0;
  for (std::size_t interface = 0; interface < pressure.size(); ++interface) {
    MPS_CHECK_NEAR(result.longwave_down_w_m2[interface],
                   source * (1.0 - std::exp(-path_from_top)), 1e-10);
    if (interface + 1 < pressure.size())
      path_from_top +=
          p.longwave_diffusivity_factor * result.optical_depth.longwave[interface];
  }
  mps::Real path_from_surface = 0.0;
  for (std::size_t reverse = pressure.size(); reverse > 0; --reverse) {
    const std::size_t interface = reverse - 1;
    MPS_CHECK_NEAR(result.longwave_up_w_m2[interface],
                   source + (boundary - source) * std::exp(-path_from_surface), 1e-10);
    if (interface > 0)
      path_from_surface +=
          p.longwave_diffusivity_factor * result.optical_depth.longwave[interface - 1];
  }
}

MPS_TEST_CASE("longwave surface reflects its unabsorbed downward flux") {
  auto p = parameters();
  p.shortwave_absorption_m2_kg = 0.0;
  p.longwave_absorption_ref_m2_kg = 5e-4;
  const std::vector<mps::Real> pressure{0.0, 100000.0};
  const std::vector<mps::Real> temperature{280.0};
  auto input = column_input(pressure, temperature, p);
  input.stellar_flux_w_m2 = 0.0;
  input.surface_emissivity = 0.3;
  const auto result = mps::gray_radiation_column(input);
  const mps::Real emitted = input.surface_emissivity * mps::kStefanBoltzmannWm2K4 *
                            std::pow(input.surface_temperature_k, 4.0);
  MPS_CHECK_NEAR(
      result.longwave_up_w_m2.back(),
      emitted + (1.0 - input.surface_emissivity) * result.longwave_down_w_m2.back(),
      1e-12);
}

MPS_TEST_CASE("opaque longwave layers approach their local source") {
  auto p = parameters();
  p.shortwave_absorption_m2_kg = 0.0;
  p.longwave_absorption_ref_m2_kg = 1e6;
  const std::vector<mps::Real> pressure{0.0, 50000.0, 100000.0};
  const std::vector<mps::Real> temperature{200.0, 260.0};
  auto input = column_input(pressure, temperature, p);
  input.stellar_flux_w_m2 = 0.0;
  const auto result = mps::gray_radiation_column(input);
  MPS_CHECK_NEAR(result.longwave_down_w_m2[1],
                 mps::kStefanBoltzmannWm2K4 * std::pow(200.0, 4.0), 1e-12);
  MPS_CHECK_NEAR(result.longwave_down_w_m2[2],
                 mps::kStefanBoltzmannWm2K4 * std::pow(260.0, 4.0), 1e-12);
  MPS_CHECK_NEAR(result.longwave_up_w_m2[0],
                 mps::kStefanBoltzmannWm2K4 * std::pow(200.0, 4.0), 1e-12);
  MPS_CHECK_NEAR(result.longwave_up_w_m2[1],
                 mps::kStefanBoltzmannWm2K4 * std::pow(260.0, 4.0), 1e-12);
}

MPS_TEST_CASE("gray column flux convergence telescopes exactly") {
  auto p = parameters();
  p.shortwave_absorption_m2_kg = 0.0;
  const std::vector<mps::Real> pressure{1000.0, 12000.0, 50000.0, 100000.0};
  const std::vector<mps::Real> temperature{210.0, 245.0, 285.0};
  const auto result =
      mps::gray_radiation_column(column_input(pressure, temperature, p));
  mps::Real sum = 0.0;
  for (const mps::Real convergence : result.radiative_convergence_w_m2)
    sum += convergence;
  MPS_CHECK_NEAR(sum, result.net_flux_w_m2.back() - result.net_flux_w_m2.front(),
                 1e-12);
}

MPS_TEST_CASE("shortwave sweeps match Beer Lambert attenuation") {
  auto p = parameters();
  p.shortwave_absorption_m2_kg = 1.5e-4;
  p.longwave_absorption_ref_m2_kg = 0.0;
  const std::vector<mps::Real> pressure{1000.0, 20000.0, 70000.0, 100000.0};
  const std::vector<mps::Real> temperature{220.0, 250.0, 280.0};
  auto input = column_input(pressure, temperature, p);
  input.cosine_solar_zenith = 0.4;
  input.surface_albedo = 0.6;
  const auto result = mps::gray_radiation_column(input);

  mps::Real direct_path = 0.0;
  for (std::size_t interface = 0; interface < pressure.size(); ++interface) {
    MPS_CHECK_NEAR(result.shortwave_down_w_m2[interface],
                   input.stellar_flux_w_m2 * input.cosine_solar_zenith *
                       std::exp(-direct_path / input.cosine_solar_zenith),
                   1e-12);
    if (interface + 1 < pressure.size())
      direct_path += result.optical_depth.shortwave[interface];
  }
  mps::Real diffuse_path = 0.0;
  const mps::Real reflected = input.surface_albedo * result.shortwave_down_w_m2.back();
  for (std::size_t reverse = pressure.size(); reverse > 0; --reverse) {
    const std::size_t interface = reverse - 1;
    MPS_CHECK_NEAR(result.shortwave_up_w_m2[interface],
                   reflected * std::exp(-p.shortwave_diffuse_factor * diffuse_path),
                   1e-12);
    if (interface > 0) diffuse_path += result.optical_depth.shortwave[interface - 1];
  }
}

MPS_TEST_CASE("shortwave input partitions into reflection and absorption") {
  auto p = parameters();
  p.shortwave_absorption_m2_kg = 9e-5;
  const std::vector<mps::Real> pressure{0.0, 10000.0, 45000.0, 100000.0};
  const std::vector<mps::Real> temperature{210.0, 250.0, 290.0};
  auto input = column_input(pressure, temperature, p);
  input.cosine_solar_zenith = 0.7;
  input.surface_albedo = 0.35;
  const auto result = mps::gray_radiation_column(input);
  mps::Real atmospheric_absorption = 0.0;
  for (const mps::Real convergence : result.shortwave_convergence_w_m2)
    atmospheric_absorption += convergence;
  const mps::Real surface_absorption =
      (1.0 - input.surface_albedo) * result.shortwave_down_w_m2.back();
  MPS_CHECK_NEAR(
      result.shortwave_down_w_m2.front(),
      result.shortwave_up_w_m2.front() + atmospheric_absorption + surface_absorption,
      1e-12);
}

MPS_TEST_CASE("night and exact terminator have no shortwave flux") {
  auto p = parameters();
  p.shortwave_absorption_m2_kg = 1e8;
  const std::vector<mps::Real> pressure{0.0, 100000.0};
  const std::vector<mps::Real> temperature{250.0};
  auto input = column_input(pressure, temperature, p);
  input.cosine_solar_zenith = 0.0;
  const auto result = mps::gray_radiation_column(input);
  for (const mps::Real flux : result.shortwave_down_w_m2)
    MPS_CHECK_NEAR(flux, 0.0, 0.0);
  for (const mps::Real flux : result.shortwave_up_w_m2) MPS_CHECK_NEAR(flux, 0.0, 0.0);
  MPS_CHECK_NEAR(result.shortwave_convergence_w_m2.front(), 0.0, 0.0);
}

MPS_TEST_CASE("small positive zenith cosine is not replaced by a floor") {
  auto p = parameters();
  p.shortwave_absorption_m2_kg = 1e-10;
  const std::vector<mps::Real> pressure{0.0, 100000.0};
  const std::vector<mps::Real> temperature{250.0};
  auto input = column_input(pressure, temperature, p);
  input.cosine_solar_zenith = 1e-12;
  const auto result = mps::gray_radiation_column(input);
  MPS_CHECK(result.shortwave_down_w_m2.front() > 0.0);
  MPS_CHECK_NEAR(result.shortwave_down_w_m2.back(), 0.0, 0.0);
}

int main() { return mps::test::run_all(); }
