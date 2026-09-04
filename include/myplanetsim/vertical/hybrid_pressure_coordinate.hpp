#pragma once

#include <span>
#include <vector>

#include "myplanetsim/core/types.hpp"

namespace mps {

struct HybridPressureCoefficients {
  std::vector<Real> a_half_pa;
  std::vector<Real> b_half;

  void validate(Real minimum_surface_pressure_pa, Real maximum_surface_pressure_pa,
                Real minimum_pressure_thickness_pa) const;
};

// The uniform sigma ramp of ADR 0007: the only coordinate family an interactive control
// request may re-resolve. It preserves the model top the preset declared and satisfies
// the ADR 0005 endpoints by construction.
[[nodiscard]] HybridPressureCoefficients uniform_sigma_coefficients(
    Real top_pressure_pa, Index levels);

// True when the given coefficients are the uniform sigma ramp at their own level count,
// so refining them stays inside the family the preset already chose.
[[nodiscard]] bool is_uniform_sigma(std::span<const Real> a_half_pa,
                                    std::span<const Real> b_half);

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
  void geometry(Real surface_pressure_pa, Real gravity_m_s2, Real gas_constant_j_kg_k,
                Real heat_capacity_cp_j_kg_k, Real reference_pressure_pa,
                HybridPressureGeometry& result) const;

 private:
  HybridPressureCoefficients coefficients_;
  Real minimum_surface_pressure_pa_;
  Real maximum_surface_pressure_pa_;
  Real minimum_pressure_thickness_pa_;
};

}  // namespace mps
