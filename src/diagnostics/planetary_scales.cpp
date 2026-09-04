#include "myplanetsim/diagnostics/planetary_scales.hpp"

#include <cmath>
#include <iomanip>
#include <locale>
#include <numbers>
#include <ostream>
#include <stdexcept>

#include "myplanetsim/core/validation.hpp"

namespace mps {
namespace {
constexpr Real kStefanBoltzmannWm2K4 = 5.670374419e-8;
}

void PlanetaryScaleInputs::validate() const {
  require_non_negative(characteristic_wind_m_s, "scales.characteristic_wind_m_s");
  require_positive(characteristic_length_m, "scales.characteristic_length_m");
  require_non_negative(buoyancy_frequency_s_1, "scales.buoyancy_frequency_s_1");
  require_positive(radiative_time_s, "scales.radiative_time_s");
  require_positive(reference_temperature_k, "scales.reference_temperature_k");
  require_non_negative(temperature_contrast_k, "scales.temperature_contrast_k");
}

PlanetaryScales compute_planetary_scales(const PlanetParameters& planet,
                                         const PlanetaryScaleInputs& inputs) {
  planet.validate();
  inputs.validate();
  PlanetaryScales result{};
  const Real absolute_rotation = std::abs(planet.rotation_rate_rad_s);
  if (absolute_rotation == 0.0) {
    return result;
  }
  const Real coriolis_reference = std::sqrt(2.0) * absolute_rotation;
  const Real scale_height =
      planet.gas_constant_j_kg_k * inputs.reference_temperature_k / planet.gravity_m_s2;
  const Real deformation_radius =
      inputs.buoyancy_frequency_s_1 * scale_height / coriolis_reference;
  result.rossby_number = {
      .value = inputs.characteristic_wind_m_s /
               (coriolis_reference * inputs.characteristic_length_m),
      .defined = true};
  result.deformation_radius_m = {.value = deformation_radius, .defined = true};
  result.burger_number = {
      .value = std::pow(deformation_radius / inputs.characteristic_length_m, 2),
      .defined = true};
  result.thermal_rossby_number = {
      .value =
          planet.gas_constant_j_kg_k * inputs.temperature_contrast_k /
          (absolute_rotation * absolute_rotation * planet.radius_m * planet.radius_m),
      .defined = true};
  result.radiative_over_rotation = {.value = inputs.radiative_time_s *
                                             absolute_rotation /
                                             (2.0 * std::numbers::pi_v<Real>),
                                    .defined = true};
  return result;
}

DefinedQuantity surface_radiative_time_scale(const Real heat_capacity_j_m2_k,
                                             const Real emissivity,
                                             const Real reference_temperature_k) {
  require_positive(heat_capacity_j_m2_k, "surface heat capacity");
  require_finite(emissivity, "surface emissivity");
  if (emissivity < 0.0 || emissivity > 1.0) {
    throw std::invalid_argument("surface emissivity must be in [0, 1]");
  }
  require_positive(reference_temperature_k, "reference temperature");
  if (emissivity == 0.0) {
    return {};
  }
  return {.value = heat_capacity_j_m2_k / (4.0 * emissivity * kStefanBoltzmannWm2K4 *
                                           std::pow(reference_temperature_k, 3)),
          .defined = true};
}

void write_defined_quantity(std::ostream& output, const std::string_view name,
                            const DefinedQuantity quantity) {
  output.imbue(std::locale::classic());
  output << std::setprecision(17) << name << ".defined = " << std::boolalpha
         << quantity.defined << '\n';
  if (quantity.defined) {
    output << name << ".value = " << quantity.value << '\n';
  }
  if (!output) {
    throw std::runtime_error("failed while writing planetary scale metadata");
  }
}

}  // namespace mps
