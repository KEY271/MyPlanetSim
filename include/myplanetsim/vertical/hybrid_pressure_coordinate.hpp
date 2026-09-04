#pragma once

#include <vector>

#include "myplanetsim/core/types.hpp"

namespace mps {

struct HybridPressureCoefficients {
  std::vector<Real> a_half_pa;
  std::vector<Real> b_half;

  void validate(Real minimum_surface_pressure_pa, Real maximum_surface_pressure_pa,
                Real minimum_pressure_thickness_pa) const;
};

struct HybridPressureGeometry {
  std::vector<Real> pressure_half_pa;
  std::vector<Real> exner_half;
  std::vector<Real> pressure_full_pa;
  std::vector<Real> exner_full;
  std::vector<Real> delta_pressure_pa;
  std::vector<Real> air_mass_kg_m2;
};

class AtmosphericHybridCoordinate {
 public:
  AtmosphericHybridCoordinate(HybridPressureCoefficients coefficients,
                              Real minimum_surface_pressure_pa,
                              Real maximum_surface_pressure_pa,
                              Real minimum_pressure_thickness_pa);

  [[nodiscard]] std::size_t levels() const noexcept;
  [[nodiscard]] Real top_pressure_pa() const noexcept;
  [[nodiscard]] const HybridPressureCoefficients& coefficients() const noexcept;
  [[nodiscard]] HybridPressureGeometry geometry(Real surface_pressure_pa,
                                                Real gravity_m_s2,
                                                Real gas_constant_j_kg_k,
                                                Real heat_capacity_cp_j_kg_k,
                                                Real reference_pressure_pa) const;

 private:
  HybridPressureCoefficients coefficients_;
  Real minimum_surface_pressure_pa_;
  Real maximum_surface_pressure_pa_;
  Real minimum_pressure_thickness_pa_;
};

}  // namespace mps
