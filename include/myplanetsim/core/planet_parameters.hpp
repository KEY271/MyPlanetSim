#pragma once

#include <optional>

#include "myplanetsim/core/types.hpp"

namespace mps {

struct PlanetParameters {
  Real radius_m;
  Real rotation_rate_rad_s;
  Real gravity_m_s2;
  Real gas_constant_j_kg_k;
  Real heat_capacity_cp_j_kg_k;
  Real reference_pressure_pa;

  [[nodiscard]] static PlanetParameters earth_like();

  void validate() const;

  [[nodiscard]] Real heat_capacity_cv_j_kg_k() const;
  [[nodiscard]] Real kappa() const;
  [[nodiscard]] std::optional<Real> rotation_period_s() const;
};

}  // namespace mps
