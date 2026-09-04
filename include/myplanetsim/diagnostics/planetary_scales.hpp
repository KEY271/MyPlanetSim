#pragma once

#include <iosfwd>
#include <string_view>

#include "myplanetsim/core/planet_parameters.hpp"
#include "myplanetsim/core/types.hpp"

namespace mps {

struct PlanetaryScaleInputs {
  Real characteristic_wind_m_s;
  Real characteristic_length_m;
  Real buoyancy_frequency_s_1;
  Real radiative_time_s;
  Real reference_temperature_k;
  Real temperature_contrast_k;

  void validate() const;
};

struct DefinedQuantity {
  Real value = 0.0;
  bool defined = false;
};

struct PlanetaryScales {
  DefinedQuantity rossby_number;
  DefinedQuantity deformation_radius_m;
  DefinedQuantity burger_number;
  DefinedQuantity thermal_rossby_number;
  DefinedQuantity radiative_over_rotation;
};

[[nodiscard]] PlanetaryScales compute_planetary_scales(
    const PlanetParameters& planet, const PlanetaryScaleInputs& inputs);
[[nodiscard]] DefinedQuantity surface_radiative_time_scale(
    Real heat_capacity_j_m2_k, Real emissivity, Real reference_temperature_k);
void write_defined_quantity(std::ostream& output, std::string_view name,
                            DefinedQuantity quantity);

}  // namespace mps
