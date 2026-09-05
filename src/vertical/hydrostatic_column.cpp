#include "myplanetsim/vertical/hydrostatic_column.hpp"

#include <cmath>
#include <stdexcept>

#include "myplanetsim/core/validation.hpp"

namespace mps {

HydrostaticColumn integrate_hydrostatic_column(
    const HybridPressureGeometry& geometry,
    const std::vector<Real>& potential_temperature_k,
    const Real heat_capacity_cp_j_kg_k, const Real gravity_m_s2,
    const Real surface_geopotential_m2_s2) {
  HydrostaticColumn result;
  integrate_hydrostatic_column(geometry, potential_temperature_k,
                               heat_capacity_cp_j_kg_k, gravity_m_s2,
                               surface_geopotential_m2_s2, result);
  return result;
}

void integrate_hydrostatic_column(const HybridPressureGeometry& geometry,
                                  const std::span<const Real> potential_temperature_k,
                                  const Real heat_capacity_cp_j_kg_k,
                                  const Real gravity_m_s2,
                                  const Real surface_geopotential_m2_s2,
                                  HydrostaticColumn& result) {
  require_positive(heat_capacity_cp_j_kg_k, "heat capacity");
  require_positive(gravity_m_s2, "gravity");
  require_finite(surface_geopotential_m2_s2, "surface geopotential");
  const std::size_t nz = geometry.exner_full.size();
  if (potential_temperature_k.size() != nz || geometry.exner_half.size() != nz + 1) {
    throw std::invalid_argument("hydrostatic geometry and theta shapes differ");
  }
  result.geopotential_half_m2_s2.resize(nz + 1);
  result.geopotential_full_m2_s2.resize(nz);
  result.height_half_m.resize(nz + 1);
  result.height_full_m.resize(nz);
  result.residual_m2_s2.resize(nz);
  result.geopotential_half_m2_s2[nz] = surface_geopotential_m2_s2;
  for (std::size_t reverse = nz; reverse > 0; --reverse) {
    const std::size_t k = reverse - 1;
    require_positive(potential_temperature_k[k], "potential temperature");
    const Real layer_integral = heat_capacity_cp_j_kg_k * potential_temperature_k[k] *
                                (geometry.exner_half[k + 1] - geometry.exner_half[k]);
    result.geopotential_half_m2_s2[k] =
        result.geopotential_half_m2_s2[k + 1] + layer_integral;
    result.geopotential_full_m2_s2[k] =
        result.geopotential_half_m2_s2[k + 1] +
        heat_capacity_cp_j_kg_k * potential_temperature_k[k] *
            (geometry.exner_half[k + 1] - geometry.exner_full[k]);
  }
  result.residual_m2_s2 = diagnose_hydrostatic_residuals(
      geometry, potential_temperature_k, heat_capacity_cp_j_kg_k, result);
  for (std::size_t k = 0; k <= nz; ++k) {
    result.height_half_m[k] = result.geopotential_half_m2_s2[k] / gravity_m_s2;
    if (k < nz) {
      result.height_full_m[k] = result.geopotential_full_m2_s2[k] / gravity_m_s2;
    }
  }
}

std::vector<Real> diagnose_hydrostatic_residuals(
    const HybridPressureGeometry& geometry,
    const std::span<const Real> potential_temperature_k,
    const Real heat_capacity_cp_j_kg_k, const HydrostaticColumn& column) {
  require_positive(heat_capacity_cp_j_kg_k, "heat capacity");
  const std::size_t nz = potential_temperature_k.size();
  if (geometry.exner_full.size() != nz || geometry.exner_half.size() != nz + 1 ||
      column.geopotential_full_m2_s2.size() != nz ||
      column.geopotential_half_m2_s2.size() != nz + 1) {
    throw std::invalid_argument("hydrostatic residual shapes differ");
  }
  std::vector<Real> residual(nz);
  for (std::size_t k = 0; k < nz; ++k) {
    require_positive(potential_temperature_k[k], "potential temperature");
    const Real reconstructed_top =
        column.geopotential_full_m2_s2[k] +
        heat_capacity_cp_j_kg_k * potential_temperature_k[k] *
            (geometry.exner_full[k] - geometry.exner_half[k]);
    residual[k] = column.geopotential_half_m2_s2[k] - reconstructed_top;
  }
  return residual;
}

}  // namespace mps
