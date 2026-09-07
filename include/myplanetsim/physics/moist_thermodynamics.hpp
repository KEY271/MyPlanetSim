#pragma once

#include "myplanetsim/core/types.hpp"

namespace mps {

struct DiluteMoistThermodynamics {
  Real gas_constant_dry_air_j_kg_k = 287.0;
  Real heat_capacity_cp_j_kg_k = 1004.0;
  Real gas_constant_water_vapor_j_kg_k = 461.5;
  Real latent_heat_vaporization_j_kg = 2.5e6;
  Real triple_point_temperature_k = 273.16;
  Real triple_point_vapor_pressure_pa = 611.657;
  Real maximum_vapor_mixing_ratio = 0.1;
  Real maximum_vapor_pressure_fraction = 0.1;

  void validate() const;
};

[[nodiscard]] Real saturation_vapor_pressure_pa(
    Real temperature_k, const DiluteMoistThermodynamics& thermodynamics);
[[nodiscard]] Real saturation_mixing_ratio(
    Real temperature_k, Real pressure_pa,
    const DiluteMoistThermodynamics& thermodynamics);
[[nodiscard]] Real saturation_mixing_ratio_temperature_derivative_k_1(
    Real temperature_k, Real pressure_pa,
    const DiluteMoistThermodynamics& thermodynamics);
[[nodiscard]] Real relative_humidity(Real vapor_mixing_ratio, Real temperature_k,
                                     Real pressure_pa,
                                     const DiluteMoistThermodynamics& thermodynamics);
void validate_dilute_moist_state(Real temperature_k, Real pressure_pa,
                                 Real vapor_mixing_ratio,
                                 const DiluteMoistThermodynamics& thermodynamics);

}  // namespace mps
