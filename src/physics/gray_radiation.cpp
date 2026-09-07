#include "myplanetsim/physics/gray_radiation.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "myplanetsim/core/validation.hpp"
#include "myplanetsim/physics/surface_energy_balance.hpp"

namespace mps {
namespace {

void validate_radiation_parameters(const RadiationParameters& parameters) {
  require_non_negative(parameters.shortwave_absorption_m2_kg,
                       "radiation.shortwave_absorption_m2_kg");
  require_non_negative(parameters.longwave_absorption_ref_m2_kg,
                       "radiation.longwave_absorption_ref_m2_kg");
  require_positive(parameters.reference_pressure_pa, "radiation.reference_pressure_pa");
  require_finite(parameters.longwave_pressure_exponent,
                 "radiation.longwave_pressure_exponent");
  if (parameters.longwave_pressure_exponent < 1.0)
    throw std::invalid_argument(
        "radiation.longwave_pressure_exponent must be at least one");
  require_positive(parameters.longwave_diffusivity_factor,
                   "radiation.longwave_diffusivity_factor");
  require_positive(parameters.shortwave_diffuse_factor,
                   "radiation.shortwave_diffuse_factor");
  require_finite(parameters.cfl, "radiation.cfl");
  if (!(parameters.cfl > 0.0 && parameters.cfl <= 1.0))
    throw std::invalid_argument("radiation.cfl must be in (0, 1]");
}

[[nodiscard]] Real pressure_power_difference(const Real top_pressure_pa,
                                             const Real bottom_pressure_pa,
                                             const Real reference_pressure_pa,
                                             const Real exponent) {
  if (top_pressure_pa == 0.0)
    return std::pow(bottom_pressure_pa / reference_pressure_pa, exponent);
  const Real top_ratio = top_pressure_pa / reference_pressure_pa;
  const Real logarithmic_ratio =
      std::log1p((bottom_pressure_pa - top_pressure_pa) / top_pressure_pa);
  return std::pow(top_ratio, exponent) * std::expm1(exponent * logarithmic_ratio);
}

}  // namespace

GrayRadiationOpticalDepth gray_radiation_optical_depth(
    const std::span<const Real> pressure_half_pa, const Real gravity_m_s2,
    const RadiationParameters& parameters) {
  GrayRadiationOpticalDepth result;
  gray_radiation_optical_depth(pressure_half_pa, gravity_m_s2, parameters, result);
  return result;
}

void gray_radiation_optical_depth(const std::span<const Real> pressure_half_pa,
                                  const Real gravity_m_s2,
                                  const RadiationParameters& parameters,
                                  GrayRadiationOpticalDepth& result) {
  validate_radiation_parameters(parameters);
  require_positive(gravity_m_s2, "gravity");
  if (pressure_half_pa.size() < 2)
    throw std::invalid_argument(
        "gray radiation requires at least two half-level pressures");

  result.shortwave.resize(pressure_half_pa.size() - 1);
  result.longwave.resize(pressure_half_pa.size() - 1);
  for (std::size_t level = 0; level + 1 < pressure_half_pa.size(); ++level) {
    const Real top = pressure_half_pa[level];
    const Real bottom = pressure_half_pa[level + 1];
    if (!std::isfinite(top) || !std::isfinite(bottom) || top < 0.0 || !(bottom > top))
      throw std::invalid_argument(
          "gray radiation half-level pressures must be finite and increasing");
    const Real pressure_thickness = bottom - top;
    result.shortwave[level] =
        parameters.shortwave_absorption_m2_kg * pressure_thickness / gravity_m_s2;
    const Real power_difference =
        pressure_power_difference(top, bottom, parameters.reference_pressure_pa,
                                  parameters.longwave_pressure_exponent);
    result.longwave[level] =
        parameters.longwave_absorption_ref_m2_kg * parameters.reference_pressure_pa /
        (gravity_m_s2 * parameters.longwave_pressure_exponent) * power_difference;
    if (!std::isfinite(result.shortwave[level]) ||
        !std::isfinite(result.longwave[level]))
      throw std::invalid_argument("gray radiation optical depth is not finite");
  }
}

GrayRadiationColumn gray_radiation_column(const GrayRadiationColumnInput& input) {
  GrayRadiationColumn result;
  gray_radiation_column(input, result);
  return result;
}

void gray_radiation_column(const GrayRadiationColumnInput& input,
                           GrayRadiationColumn& result) {
  const std::size_t levels = input.temperature_k.size();
  if (levels == 0 || input.pressure_half_pa.size() != levels + 1)
    throw std::invalid_argument("gray radiation column shape mismatch");
  require_positive(input.surface_temperature_k, "surface temperature");
  require_positive(input.gravity_m_s2, "gravity");
  require_non_negative(input.stellar_flux_w_m2, "stellar flux");
  require_finite(input.cosine_solar_zenith, "cosine solar zenith");
  if (input.cosine_solar_zenith < 0.0 || input.cosine_solar_zenith > 1.0)
    throw std::invalid_argument("cosine solar zenith must be in [0, 1]");
  require_finite(input.surface_albedo, "surface albedo");
  require_finite(input.surface_emissivity, "surface emissivity");
  if (input.surface_albedo < 0.0 || input.surface_albedo > 1.0 ||
      input.surface_emissivity < 0.0 || input.surface_emissivity > 1.0)
    throw std::invalid_argument("surface albedo and emissivity must be in [0, 1]");
  for (const Real temperature : input.temperature_k)
    require_positive(temperature, "atmospheric temperature");

  gray_radiation_optical_depth(input.pressure_half_pa, input.gravity_m_s2,
                               input.parameters, result.optical_depth);
  result.shortwave_down_w_m2.resize(levels + 1);
  result.shortwave_up_w_m2.resize(levels + 1);
  result.longwave_down_w_m2.resize(levels + 1);
  result.longwave_up_w_m2.resize(levels + 1);
  result.net_flux_w_m2.resize(levels + 1);
  result.shortwave_convergence_w_m2.resize(levels);
  result.longwave_convergence_w_m2.resize(levels);
  result.radiative_convergence_w_m2.resize(levels);

  const Real incoming_shortwave = input.stellar_flux_w_m2 * input.cosine_solar_zenith;
  result.shortwave_down_w_m2[0] = incoming_shortwave;
  if (input.cosine_solar_zenith == 0.0) {
    std::fill(result.shortwave_down_w_m2.begin(), result.shortwave_down_w_m2.end(),
              0.0);
    std::fill(result.shortwave_up_w_m2.begin(), result.shortwave_up_w_m2.end(), 0.0);
  } else {
    for (std::size_t level = 0; level < levels; ++level) {
      result.shortwave_down_w_m2[level + 1] =
          result.shortwave_down_w_m2[level] *
          std::exp(-result.optical_depth.shortwave[level] / input.cosine_solar_zenith);
    }
    result.shortwave_up_w_m2[levels] =
        input.surface_albedo * result.shortwave_down_w_m2[levels];
    for (std::size_t reverse = levels; reverse > 0; --reverse) {
      const std::size_t level = reverse - 1;
      result.shortwave_up_w_m2[level] =
          result.shortwave_up_w_m2[level + 1] *
          std::exp(-input.parameters.shortwave_diffuse_factor *
                   result.optical_depth.shortwave[level]);
    }
  }

  result.longwave_down_w_m2[0] = 0.0;
  for (std::size_t level = 0; level < levels; ++level) {
    const Real path = input.parameters.longwave_diffusivity_factor *
                      result.optical_depth.longwave[level];
    const Real transmission = std::exp(-path);
    const Real absorptivity = -std::expm1(-path);
    const Real temperature_squared =
        input.temperature_k[level] * input.temperature_k[level];
    const Real source =
        kStefanBoltzmannWm2K4 * temperature_squared * temperature_squared;
    result.longwave_down_w_m2[level + 1] =
        transmission * result.longwave_down_w_m2[level] + absorptivity * source;
  }

  const Real surface_temperature_squared =
      input.surface_temperature_k * input.surface_temperature_k;
  result.longwave_up_w_m2[levels] =
      input.surface_emissivity * kStefanBoltzmannWm2K4 * surface_temperature_squared *
          surface_temperature_squared +
      (1.0 - input.surface_emissivity) * result.longwave_down_w_m2[levels];
  for (std::size_t reverse = levels; reverse > 0; --reverse) {
    const std::size_t level = reverse - 1;
    const Real path = input.parameters.longwave_diffusivity_factor *
                      result.optical_depth.longwave[level];
    const Real transmission = std::exp(-path);
    const Real absorptivity = -std::expm1(-path);
    const Real temperature_squared =
        input.temperature_k[level] * input.temperature_k[level];
    const Real source =
        kStefanBoltzmannWm2K4 * temperature_squared * temperature_squared;
    result.longwave_up_w_m2[level] =
        transmission * result.longwave_up_w_m2[level + 1] + absorptivity * source;
  }

  for (std::size_t interface = 0; interface <= levels; ++interface) {
    result.net_flux_w_m2[interface] =
        result.longwave_up_w_m2[interface] - result.longwave_down_w_m2[interface] +
        result.shortwave_up_w_m2[interface] - result.shortwave_down_w_m2[interface];
  }
  for (std::size_t level = 0; level < levels; ++level) {
    result.shortwave_convergence_w_m2[level] =
        (result.shortwave_up_w_m2[level + 1] - result.shortwave_down_w_m2[level + 1]) -
        (result.shortwave_up_w_m2[level] - result.shortwave_down_w_m2[level]);
    result.longwave_convergence_w_m2[level] =
        (result.longwave_up_w_m2[level + 1] - result.longwave_down_w_m2[level + 1]) -
        (result.longwave_up_w_m2[level] - result.longwave_down_w_m2[level]);
    result.radiative_convergence_w_m2[level] =
        result.shortwave_convergence_w_m2[level] +
        result.longwave_convergence_w_m2[level];
  }
}

GrayRadiativeColumnTendency gray_radiative_column_tendency(
    const GrayRadiativeColumnSourceInput& input) {
  GrayRadiativeColumnTendency result;
  gray_radiative_column_tendency(input, result);
  return result;
}

void gray_radiative_column_tendency(const GrayRadiativeColumnSourceInput& input,
                                    GrayRadiativeColumnTendency& result) {
  const std::size_t levels = input.radiation.temperature_k.size();
  if (levels == 0 || input.exner_full.size() != levels)
    throw std::invalid_argument("gray radiative source shape mismatch");
  require_positive(input.heat_capacity_cp_j_kg_k, "heat capacity cp");
  require_positive(input.surface_heat_capacity_j_m2_k, "surface heat capacity");
  require_non_negative(input.air_exchange_coefficient_w_m2_k,
                       "air exchange coefficient");
  require_finite(input.internal_heat_flux_w_m2, "internal heat flux");
  for (const Real exner : input.exner_full) require_positive(exner, "full-level Exner");

  gray_radiation_column(input.radiation, result.fluxes);
  result.potential_temperature_mass_k_kg_m2_s.resize(levels);
  for (std::size_t level = 0; level < levels; ++level) {
    result.potential_temperature_mass_k_kg_m2_s[level] =
        result.fluxes.radiative_convergence_w_m2[level] /
        (input.heat_capacity_cp_j_kg_k * input.exner_full[level]);
  }

  const Real sensible =
      input.air_exchange_coefficient_w_m2_k *
      (input.radiation.surface_temperature_k - input.radiation.temperature_k.back());
  result.sensible_to_atmosphere_w_m2 = sensible;
  result.potential_temperature_mass_k_kg_m2_s.back() +=
      sensible / (input.heat_capacity_cp_j_kg_k * input.exner_full.back());
  result.surface_storage_rate_w_m2 =
      -result.fluxes.net_flux_w_m2.back() + input.internal_heat_flux_w_m2 - sensible;
  result.surface_temperature_k_s =
      result.surface_storage_rate_w_m2 / input.surface_heat_capacity_j_m2_k;

  // Propagate L1 norms of temperature sensitivities through the positive-coefficient
  // longwave sweeps. Triangle inequalities at net-flux differences make this a
  // conservative O(K) bound on every temperature-tendency Jacobian row.
  std::vector<Real> downward_sensitivity(levels + 1, 0.0);
  std::vector<Real> upward_sensitivity(levels + 1, 0.0);
  for (std::size_t level = 0; level < levels; ++level) {
    const Real path = input.radiation.parameters.longwave_diffusivity_factor *
                      result.fluxes.optical_depth.longwave[level];
    const Real transmission = std::exp(-path);
    const Real absorptivity = -std::expm1(-path);
    const Real temperature = input.radiation.temperature_k[level];
    const Real source_derivative =
        4.0 * kStefanBoltzmannWm2K4 * temperature * temperature * temperature;
    downward_sensitivity[level + 1] =
        transmission * downward_sensitivity[level] + absorptivity * source_derivative;
  }
  const Real surface_temperature = input.radiation.surface_temperature_k;
  const Real surface_source_derivative = 4.0 * kStefanBoltzmannWm2K4 *
                                         surface_temperature * surface_temperature *
                                         surface_temperature;
  upward_sensitivity[levels] =
      input.radiation.surface_emissivity * surface_source_derivative +
      (1.0 - input.radiation.surface_emissivity) * downward_sensitivity[levels];
  for (std::size_t reverse = levels; reverse > 0; --reverse) {
    const std::size_t level = reverse - 1;
    const Real path = input.radiation.parameters.longwave_diffusivity_factor *
                      result.fluxes.optical_depth.longwave[level];
    const Real transmission = std::exp(-path);
    const Real absorptivity = -std::expm1(-path);
    const Real temperature = input.radiation.temperature_k[level];
    const Real source_derivative =
        4.0 * kStefanBoltzmannWm2K4 * temperature * temperature * temperature;
    upward_sensitivity[level] =
        transmission * upward_sensitivity[level + 1] + absorptivity * source_derivative;
  }

  Real rate_bound = 0.0;
  for (std::size_t level = 0; level < levels; ++level) {
    const Real mass = (input.radiation.pressure_half_pa[level + 1] -
                       input.radiation.pressure_half_pa[level]) /
                      input.radiation.gravity_m_s2;
    const Real flux_sensitivity =
        upward_sensitivity[level] + downward_sensitivity[level] +
        upward_sensitivity[level + 1] + downward_sensitivity[level + 1];
    Real row_bound = flux_sensitivity / (input.heat_capacity_cp_j_kg_k * mass);
    if (level + 1 == levels)
      row_bound += 2.0 * input.air_exchange_coefficient_w_m2_k /
                   (input.heat_capacity_cp_j_kg_k * mass);
    rate_bound = std::max(rate_bound, row_bound);
  }
  const Real surface_row_bound =
      (upward_sensitivity[levels] + downward_sensitivity[levels] +
       2.0 * input.air_exchange_coefficient_w_m2_k) /
      input.surface_heat_capacity_j_m2_k;
  rate_bound = std::max(rate_bound, surface_row_bound);
  if (!std::isfinite(rate_bound))
    throw std::invalid_argument("gray radiation temperature rate bound is not finite");
  result.temperature_rate_bound_s_1 = rate_bound;
  result.stable_time_step_s = rate_bound == 0.0
                                  ? std::numeric_limits<Real>::infinity()
                                  : input.radiation.parameters.cfl / rate_bound;
}

}  // namespace mps
