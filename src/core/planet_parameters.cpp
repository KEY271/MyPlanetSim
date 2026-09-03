#include "myplanetsim/core/planet_parameters.hpp"

#include <cmath>
#include <numbers>
#include <stdexcept>

#include "myplanetsim/core/validation.hpp"

namespace mps {

PlanetParameters PlanetParameters::earth_like() {
  return PlanetParameters{
      .radius_m = 6.37122e6,
      .rotation_rate_rad_s = 7.292115e-5,
      .gravity_m_s2 = 9.80616,
      .gas_constant_j_kg_k = 287.0,
      .heat_capacity_cp_j_kg_k = 1004.5,
      .reference_pressure_pa = 100000.0,
  };
}

void PlanetParameters::validate() const {
  require_positive(radius_m, "planet.radius_m");
  require_finite(rotation_rate_rad_s, "planet.rotation_rate_rad_s");
  require_positive(gravity_m_s2, "planet.gravity_m_s2");
  require_positive(gas_constant_j_kg_k, "planet.gas_constant_j_kg_k");
  require_positive(heat_capacity_cp_j_kg_k, "planet.heat_capacity_cp_j_kg_k");
  require_positive(reference_pressure_pa, "planet.reference_pressure_pa");

  if (heat_capacity_cp_j_kg_k <= gas_constant_j_kg_k) {
    throw std::invalid_argument(
        "planet.heat_capacity_cp_j_kg_k must be greater than "
        "planet.gas_constant_j_kg_k");
  }
}

Real PlanetParameters::heat_capacity_cv_j_kg_k() const {
  return heat_capacity_cp_j_kg_k - gas_constant_j_kg_k;
}

Real PlanetParameters::kappa() const {
  return gas_constant_j_kg_k / heat_capacity_cp_j_kg_k;
}

std::optional<Real> PlanetParameters::rotation_period_s() const {
  if (rotation_rate_rad_s == 0.0) {
    return std::nullopt;
  }
  return 2.0 * std::numbers::pi_v<Real> / std::abs(rotation_rate_rad_s);
}

}  // namespace mps
