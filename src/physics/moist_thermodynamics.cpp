#include "myplanetsim/physics/moist_thermodynamics.hpp"

#include <cmath>
#include <stdexcept>

#include "myplanetsim/core/validation.hpp"

namespace mps {

void DiluteMoistThermodynamics::validate() const {
  require_positive(gas_constant_dry_air_j_kg_k, "moist dry-air gas constant");
  require_positive(heat_capacity_cp_j_kg_k, "moist heat capacity cp");
  if (!(heat_capacity_cp_j_kg_k > gas_constant_dry_air_j_kg_k))
    throw std::invalid_argument("moist cp must exceed the dry-air gas constant");
  require_positive(gas_constant_water_vapor_j_kg_k, "water-vapor gas constant");
  require_positive(latent_heat_vaporization_j_kg, "latent heat of vaporization");
  require_positive(triple_point_temperature_k, "triple-point temperature");
  require_positive(triple_point_vapor_pressure_pa, "triple-point vapor pressure");
  require_positive(maximum_vapor_mixing_ratio, "maximum vapor mixing ratio");
  require_positive(maximum_vapor_pressure_fraction, "maximum vapor pressure fraction");
  if (maximum_vapor_mixing_ratio > 1.0 || maximum_vapor_pressure_fraction > 1.0)
    throw std::invalid_argument("dilute-water limits must not exceed one");
}

Real saturation_vapor_pressure_pa(const Real temperature_k,
                                  const DiluteMoistThermodynamics& thermodynamics) {
  thermodynamics.validate();
  require_positive(temperature_k, "moist temperature");
  const Real exponent =
      thermodynamics.latent_heat_vaporization_j_kg /
      thermodynamics.gas_constant_water_vapor_j_kg_k *
      (1.0 / thermodynamics.triple_point_temperature_k - 1.0 / temperature_k);
  const Real result =
      thermodynamics.triple_point_vapor_pressure_pa * std::exp(exponent);
  if (!(result > 0.0) || !std::isfinite(result))
    throw std::invalid_argument("saturation vapor pressure is outside the domain");
  return result;
}

Real saturation_mixing_ratio(const Real temperature_k, const Real pressure_pa,
                             const DiluteMoistThermodynamics& thermodynamics) {
  require_positive(pressure_pa, "moist pressure");
  const Real epsilon = thermodynamics.gas_constant_dry_air_j_kg_k /
                       thermodynamics.gas_constant_water_vapor_j_kg_k;
  return epsilon * saturation_vapor_pressure_pa(temperature_k, thermodynamics) /
         pressure_pa;
}

Real saturation_mixing_ratio_temperature_derivative_k_1(
    const Real temperature_k, const Real pressure_pa,
    const DiluteMoistThermodynamics& thermodynamics) {
  const Real saturation =
      saturation_mixing_ratio(temperature_k, pressure_pa, thermodynamics);
  return saturation * thermodynamics.latent_heat_vaporization_j_kg /
         (thermodynamics.gas_constant_water_vapor_j_kg_k * temperature_k *
          temperature_k);
}

Real relative_humidity(const Real vapor_mixing_ratio, const Real temperature_k,
                       const Real pressure_pa,
                       const DiluteMoistThermodynamics& thermodynamics) {
  if (vapor_mixing_ratio < 0.0 || !std::isfinite(vapor_mixing_ratio))
    throw std::invalid_argument("vapor mixing ratio must be finite and nonnegative");
  return vapor_mixing_ratio /
         saturation_mixing_ratio(temperature_k, pressure_pa, thermodynamics);
}

void validate_dilute_moist_state(const Real temperature_k, const Real pressure_pa,
                                 const Real vapor_mixing_ratio,
                                 const DiluteMoistThermodynamics& thermodynamics) {
  thermodynamics.validate();
  require_positive(temperature_k, "moist temperature");
  require_positive(pressure_pa, "moist pressure");
  if (vapor_mixing_ratio < 0.0 || !std::isfinite(vapor_mixing_ratio))
    throw std::invalid_argument("vapor mixing ratio must be finite and nonnegative");
  const Real pressure_fraction =
      saturation_vapor_pressure_pa(temperature_k, thermodynamics) / pressure_pa;
  if (vapor_mixing_ratio > thermodynamics.maximum_vapor_mixing_ratio ||
      pressure_fraction > thermodynamics.maximum_vapor_pressure_fraction)
    throw std::domain_error("state is outside the configured dilute-water domain");
}

}  // namespace mps
