#include <algorithm>
#include <array>
#include <cmath>
#include <span>
#include <vector>

#include "myplanetsim/numerics/explicit_steppers.hpp"
#include "myplanetsim/physics/gray_radiation.hpp"
#include "myplanetsim/physics/surface_energy_balance.hpp"
#include "support/test.hpp"

namespace {

[[nodiscard]] mps::RadiationParameters radiation_parameters() {
  return mps::RadiationParameters{
      .shortwave_absorption_m2_kg = 0.0,
      .longwave_absorption_ref_m2_kg = 3e-4,
      .reference_pressure_pa = 100000.0,
      .longwave_pressure_exponent = 1.0,
      .longwave_diffusivity_factor = 1.66,
      .shortwave_diffuse_factor = 1.66,
      .cfl = 0.5,
  };
}

[[nodiscard]] mps::GrayRadiativeColumnSourceInput source_input(
    const std::vector<mps::Real>& pressure, const std::vector<mps::Real>& temperature,
    const std::vector<mps::Real>& exner, const mps::Real surface_temperature,
    const mps::RadiationParameters& parameters) {
  return mps::GrayRadiativeColumnSourceInput{
      .radiation =
          mps::GrayRadiationColumnInput{
              .pressure_half_pa = pressure,
              .temperature_k = temperature,
              .surface_temperature_k = surface_temperature,
              .gravity_m_s2 = 10.0,
              .stellar_flux_w_m2 = 0.0,
              .cosine_solar_zenith = 0.0,
              .surface_albedo = 0.0,
              .surface_emissivity = 1.0,
              .parameters = parameters,
          },
      .exner_full = exner,
      .heat_capacity_cp_j_kg_k = 1000.0,
      .surface_heat_capacity_j_m2_k = 1e7,
      .air_exchange_coefficient_w_m2_k = 0.0,
      .internal_heat_flux_w_m2 = 0.0,
  };
}

[[nodiscard]] std::vector<mps::Real> temperature_tendency(
    const std::vector<mps::Real>& pressure, const std::vector<mps::Real>& exner,
    const std::vector<mps::Real>& state, const mps::RadiationParameters& parameters,
    const mps::Real exchange) {
  const std::size_t levels = exner.size();
  std::vector<mps::Real> temperature(state.begin(), state.begin() + levels);
  auto input = source_input(pressure, temperature, exner, state.back(), parameters);
  input.air_exchange_coefficient_w_m2_k = exchange;
  const auto source = mps::gray_radiative_column_tendency(input);
  std::vector<mps::Real> result(levels + 1);
  for (std::size_t level = 0; level < levels; ++level) {
    const mps::Real mass =
        (pressure[level + 1] - pressure[level]) / input.radiation.gravity_m_s2;
    result[level] =
        source.potential_temperature_mass_k_kg_m2_s[level] * exner[level] / mass;
  }
  result.back() = source.surface_temperature_k_s;
  return result;
}

}  // namespace

MPS_TEST_CASE("radiative convergence and sensible heat use the same reservoirs") {
  auto parameters = radiation_parameters();
  parameters.longwave_absorption_ref_m2_kg = 0.0;
  const std::vector<mps::Real> pressure{0.0, 100000.0};
  const std::vector<mps::Real> temperature{280.0};
  const std::vector<mps::Real> exner{0.9};
  auto input = source_input(pressure, temperature, exner, 300.0, parameters);
  input.radiation.surface_emissivity = 0.0;
  input.air_exchange_coefficient_w_m2_k = 10.0;
  const auto source = mps::gray_radiative_column_tendency(input);
  MPS_CHECK_NEAR(source.sensible_to_atmosphere_w_m2, 200.0, 1e-13);
  MPS_CHECK_NEAR(source.potential_temperature_mass_k_kg_m2_s.front() *
                     input.heat_capacity_cp_j_kg_k * exner.front(),
                 200.0, 1e-12);
  MPS_CHECK_NEAR(source.surface_storage_rate_w_m2, -200.0, 1e-12);
  MPS_CHECK_NEAR(source.surface_temperature_k_s,
                 -200.0 / input.surface_heat_capacity_j_m2_k, 1e-20);
}

MPS_TEST_CASE("one-layer gray greenhouse equilibrium has zero source") {
  constexpr mps::Real absorption = 0.6;
  constexpr mps::Real incoming = 240.0;
  auto parameters = radiation_parameters();
  parameters.longwave_absorption_ref_m2_kg =
      -std::log(1.0 - absorption) * 10.0 /
      (parameters.longwave_diffusivity_factor * 100000.0);
  const mps::Real surface_emission = incoming / (1.0 - absorption / 2.0);
  const mps::Real surface_temperature =
      std::pow(surface_emission / mps::kStefanBoltzmannWm2K4, 0.25);
  const mps::Real air_temperature = surface_temperature / std::pow(2.0, 0.25);
  const std::vector<mps::Real> pressure{0.0, 100000.0};
  const std::vector<mps::Real> temperature{air_temperature};
  const std::vector<mps::Real> exner{1.0};
  auto input =
      source_input(pressure, temperature, exner, surface_temperature, parameters);
  input.radiation.stellar_flux_w_m2 = incoming;
  input.radiation.cosine_solar_zenith = 1.0;
  const auto source = mps::gray_radiative_column_tendency(input);
  MPS_CHECK_NEAR(source.potential_temperature_mass_k_kg_m2_s.front(), 0.0, 1e-15);
  MPS_CHECK_NEAR(source.surface_temperature_k_s, 0.0, 1e-18);
}

MPS_TEST_CASE("O K temperature rate bound encloses finite difference Jacobian") {
  auto parameters = radiation_parameters();
  parameters.shortwave_absorption_m2_kg = 8e-5;
  const std::vector<mps::Real> pressure{1000.0, 40000.0, 100000.0};
  const std::vector<mps::Real> exner{0.7, 0.95};
  const std::vector<mps::Real> state{230.0, 285.0, 305.0};
  constexpr mps::Real exchange = 12.0;
  std::array<mps::Real, 3> row_sum{};
  constexpr mps::Real epsilon = 1e-3;
  for (std::size_t variable = 0; variable < state.size(); ++variable) {
    auto upper = state;
    auto lower = state;
    upper[variable] += epsilon;
    lower[variable] -= epsilon;
    const auto upper_rate =
        temperature_tendency(pressure, exner, upper, parameters, exchange);
    const auto lower_rate =
        temperature_tendency(pressure, exner, lower, parameters, exchange);
    for (std::size_t row = 0; row < row_sum.size(); ++row)
      row_sum[row] += std::abs((upper_rate[row] - lower_rate[row]) / (2.0 * epsilon));
  }
  std::vector<mps::Real> temperature(state.begin(), state.end() - 1);
  auto input = source_input(pressure, temperature, exner, state.back(), parameters);
  input.air_exchange_coefficient_w_m2_k = exchange;
  const auto source = mps::gray_radiative_column_tendency(input);
  const mps::Real measured = *std::max_element(row_sum.begin(), row_sum.end());
  MPS_CHECK(source.temperature_rate_bound_s_1 >= measured * (1.0 - 1e-8));
  MPS_CHECK_NEAR(source.stable_time_step_s,
                 parameters.cfl / source.temperature_rate_bound_s_1, 1e-12);
}

MPS_TEST_CASE("surface-only cooling has third-order SSPRK3 convergence") {
  auto parameters = radiation_parameters();
  parameters.longwave_absorption_ref_m2_kg = 0.0;
  const std::vector<mps::Real> pressure{0.0, 100000.0};
  const std::vector<mps::Real> temperature{250.0};
  const std::vector<mps::Real> exner{1.0};
  constexpr mps::Real initial_temperature = 300.0;
  constexpr mps::Real final_time = 1e6;
  constexpr mps::Real capacity = 1e7;
  const mps::Real exact =
      std::pow(std::pow(initial_temperature, -3.0) +
                   3.0 * mps::kStefanBoltzmannWm2K4 * final_time / capacity,
               -1.0 / 3.0);
  std::array<mps::Real, 3> errors{};
  for (std::size_t refinement = 0; refinement < errors.size(); ++refinement) {
    const std::size_t steps = 10U << refinement;
    const mps::Real time_step = final_time / static_cast<mps::Real>(steps);
    std::vector<mps::Real> state{initial_temperature};
    mps::SspRk3 stepper(1);
    const auto rhs = [&](const mps::Real, const std::span<const mps::Real> value,
                         const std::span<mps::Real> rate) {
      auto input =
          source_input(pressure, temperature, exner, value.front(), parameters);
      input.surface_heat_capacity_j_m2_k = capacity;
      rate.front() = mps::gray_radiative_column_tendency(input).surface_temperature_k_s;
    };
    for (std::size_t step = 0; step < steps; ++step)
      stepper.step(static_cast<mps::Real>(step) * time_step, time_step, state, rhs);
    errors[refinement] = std::abs(state.front() - exact);
  }
  MPS_CHECK(std::log2(errors[0] / errors[1]) >= 2.7);
  MPS_CHECK(std::log2(errors[1] / errors[2]) >= 2.7);
}

MPS_TEST_CASE("radiative safety bound causes an oversized step to be retried") {
  auto parameters = radiation_parameters();
  const std::vector<mps::Real> pressure{0.0, 100000.0};
  const std::vector<mps::Real> temperature{280.0};
  const std::vector<mps::Real> exner{1.0};
  const auto input = source_input(pressure, temperature, exner, 320.0, parameters);
  const auto source = mps::gray_radiative_column_tendency(input);
  mps::Real attempted = 4.0 * source.stable_time_step_s;
  std::size_t retries = 0;
  while (attempted > source.stable_time_step_s) {
    attempted *= 0.5;
    ++retries;
  }
  MPS_CHECK_EQ(retries, 2U);
  MPS_CHECK(attempted <= source.stable_time_step_s);
}

int main() { return mps::test::run_all(); }
