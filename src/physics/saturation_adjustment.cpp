#include "myplanetsim/physics/saturation_adjustment.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace mps {
namespace {

[[nodiscard]] Real layer_enthalpy(const Real mass, const Real temperature,
                                  const Real vapor_mixing_ratio,
                                  const DiluteMoistThermodynamics& thermodynamics) {
  return mass * (thermodynamics.heat_capacity_cp_j_kg_k * temperature +
                 thermodynamics.latent_heat_vaporization_j_kg * vapor_mixing_ratio);
}

}  // namespace

SaturationAdjustmentDiagnostics adjust_saturation_column(
    const std::span<const Real> air_mass_kg_m2, const std::span<const Real> pressure_pa,
    const std::span<const Real> exner_full,
    const std::span<Real> potential_temperature_mass_k_kg_m2,
    const std::span<Real> vapor_mass_kg_m2,
    const DiluteMoistThermodynamics& thermodynamics) {
  thermodynamics.validate();
  const auto levels = air_mass_kg_m2.size();
  if (levels == 0 || pressure_pa.size() != levels || exner_full.size() != levels ||
      potential_temperature_mass_k_kg_m2.size() != levels ||
      vapor_mass_kg_m2.size() != levels)
    throw std::invalid_argument("saturation-adjustment column shape mismatch");

  SaturationAdjustmentDiagnostics result;
  Real enthalpy_before = 0.0;
  Real enthalpy_after = 0.0;
  for (std::size_t level = 0; level < levels; ++level) {
    const Real mass = air_mass_kg_m2[level];
    const Real exner = exner_full[level];
    if (!(mass > 0.0) || !std::isfinite(mass) || !(exner > 0.0) ||
        !std::isfinite(exner))
      throw std::invalid_argument("saturation-adjustment geometry is invalid");
    const Real theta = potential_temperature_mass_k_kg_m2[level] / mass;
    const Real temperature = theta * exner;
    const Real vapor = vapor_mass_kg_m2[level] / mass;
    validate_dilute_moist_state(temperature, pressure_pa[level], vapor, thermodynamics);
    enthalpy_before += layer_enthalpy(mass, temperature, vapor, thermodynamics);
    const Real initial_saturation =
        saturation_mixing_ratio(temperature, pressure_pa[level], thermodynamics);
    const Real saturation_tolerance =
        1e-14 + 1e-12 * std::max(vapor, initial_saturation);
    result.maximum_supersaturation =
        std::max(result.maximum_supersaturation, vapor - initial_saturation);
    if (vapor - initial_saturation <= saturation_tolerance) {
      enthalpy_after += layer_enthalpy(mass, temperature, vapor, thermodynamics);
      result.maximum_vapor_mixing_ratio =
          std::max(result.maximum_vapor_mixing_ratio, vapor);
      result.maximum_vapor_pressure_fraction =
          std::max(result.maximum_vapor_pressure_fraction,
                   saturation_vapor_pressure_pa(temperature, thermodynamics) /
                       pressure_pa[level]);
      continue;
    }

    const Real warming_per_mixing_ratio = thermodynamics.latent_heat_vaporization_j_kg /
                                          thermodynamics.heat_capacity_cp_j_kg_k;
    const auto residual = [&](const Real condensed) {
      return vapor - condensed -
             saturation_mixing_ratio(temperature + warming_per_mixing_ratio * condensed,
                                     pressure_pa[level], thermodynamics);
    };
    Real lower = 0.0;
    Real upper = vapor;
    Real condensed =
        std::clamp((vapor - initial_saturation) /
                       (1.0 + warming_per_mixing_ratio *
                                  saturation_mixing_ratio_temperature_derivative_k_1(
                                      temperature, pressure_pa[level], thermodynamics)),
                   lower, upper);
    for (std::size_t iteration = 0; iteration < 100; ++iteration) {
      const Real value = residual(condensed);
      ++result.nonlinear_iterations;
      const Real scale = std::max({1e-16, vapor, initial_saturation});
      if (std::abs(value) <= 1e-14 + 1e-12 * scale) break;
      if (value > 0.0)
        lower = condensed;
      else
        upper = condensed;
      const Real trial_temperature = temperature + warming_per_mixing_ratio * condensed;
      const Real derivative =
          -1.0 - warming_per_mixing_ratio *
                     saturation_mixing_ratio_temperature_derivative_k_1(
                         trial_temperature, pressure_pa[level], thermodynamics);
      const Real newton = condensed - value / derivative;
      condensed = newton > lower && newton < upper ? newton : 0.5 * (lower + upper);
      if (iteration == 99)
        throw std::runtime_error("saturation adjustment did not converge");
    }
    const Real final_temperature = temperature + warming_per_mixing_ratio * condensed;
    const Real final_vapor = vapor - condensed;
    validate_dilute_moist_state(final_temperature, pressure_pa[level], final_vapor,
                                thermodynamics);
    potential_temperature_mass_k_kg_m2[level] = mass * final_temperature / exner;
    vapor_mass_kg_m2[level] = mass * final_vapor;
    result.condensed_water_kg_m2 += mass * condensed;
    ++result.adjusted_layer_count;
    result.maximum_vapor_mixing_ratio =
        std::max(result.maximum_vapor_mixing_ratio, final_vapor);
    result.maximum_vapor_pressure_fraction =
        std::max(result.maximum_vapor_pressure_fraction,
                 saturation_vapor_pressure_pa(final_temperature, thermodynamics) /
                     pressure_pa[level]);
    enthalpy_after +=
        layer_enthalpy(mass, final_temperature, final_vapor, thermodynamics);
  }
  result.enthalpy_change_j_m2 = enthalpy_after - enthalpy_before;
  return result;
}

}  // namespace mps
