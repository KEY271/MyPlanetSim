#pragma once

#include <vector>

#include "myplanetsim/core/types.hpp"
#include "myplanetsim/vertical/hybrid_pressure_coordinate.hpp"

namespace mps {

struct HydrostaticColumn {
  std::vector<Real> geopotential_half_m2_s2;
  std::vector<Real> geopotential_full_m2_s2;
  std::vector<Real> height_half_m;
  std::vector<Real> height_full_m;
  std::vector<Real> residual_m2_s2;
};

[[nodiscard]] HydrostaticColumn integrate_hydrostatic_column(
    const HybridPressureGeometry& geometry,
    const std::vector<Real>& potential_temperature_k, Real heat_capacity_cp_j_kg_k,
    Real gravity_m_s2, Real surface_geopotential_m2_s2);

}  // namespace mps
