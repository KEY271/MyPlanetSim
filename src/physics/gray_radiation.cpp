#include "myplanetsim/physics/gray_radiation.hpp"

#include <cmath>
#include <stdexcept>

#include "myplanetsim/core/validation.hpp"

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
  validate_radiation_parameters(parameters);
  require_positive(gravity_m_s2, "gravity");
  if (pressure_half_pa.size() < 2)
    throw std::invalid_argument(
        "gray radiation requires at least two half-level pressures");

  GrayRadiationOpticalDepth result;
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
  return result;
}

}  // namespace mps
