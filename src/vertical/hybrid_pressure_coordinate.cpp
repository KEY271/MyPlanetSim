#include "myplanetsim/vertical/hybrid_pressure_coordinate.hpp"

#include <cmath>
#include <stdexcept>

#include "myplanetsim/core/validation.hpp"

namespace mps {
namespace {

[[nodiscard]] Real stable_exner_mean(const Real upper_pressure_pa,
                                     const Real lower_pressure_pa, const Real kappa,
                                     const Real reference_pressure_pa) {
  const Real log_ratio = std::log(lower_pressure_pa / upper_pressure_pa);
  const Real upper = std::pow(upper_pressure_pa / reference_pressure_pa, kappa);
  if (std::abs(log_ratio) < 1.0e-8) {
    return upper * std::expm1(kappa * log_ratio) / (kappa * log_ratio);
  }
  const Real lower = std::pow(lower_pressure_pa / reference_pressure_pa, kappa);
  return (lower - upper) / (kappa * log_ratio);
}

}  // namespace

HybridPressureCoefficients uniform_sigma_coefficients(const Real top_pressure_pa,
                                                      const Index levels) {
  require_positive(top_pressure_pa, "hybrid model top pressure");
  if (levels < 1) {
    throw std::invalid_argument("hybrid coefficient level count must be positive");
  }
  const auto count = static_cast<std::size_t>(levels) + 1;
  HybridPressureCoefficients coefficients;
  coefficients.a_half_pa.reserve(count);
  coefficients.b_half.reserve(count);
  for (std::size_t k = 0; k < count; ++k) {
    const Real sigma = static_cast<Real>(k) / static_cast<Real>(levels);
    coefficients.a_half_pa.push_back(top_pressure_pa * (1.0 - sigma));
    coefficients.b_half.push_back(sigma);
  }
  // The endpoints are the ADR 0005 contract, so they are set exactly rather than left
  // to the rounding of the ramp.
  coefficients.a_half_pa.front() = top_pressure_pa;
  coefficients.a_half_pa.back() = 0.0;
  coefficients.b_half.front() = 0.0;
  coefficients.b_half.back() = 1.0;
  return coefficients;
}

bool is_uniform_sigma(const std::span<const Real> a_half_pa,
                      const std::span<const Real> b_half) {
  if (a_half_pa.size() < 2 || a_half_pa.size() != b_half.size()) {
    return false;
  }
  const auto levels = static_cast<Index>(a_half_pa.size() - 1);
  const Real top_pressure_pa = a_half_pa.front();
  if (!(top_pressure_pa > 0.0)) {
    return false;
  }
  const auto expected = uniform_sigma_coefficients(top_pressure_pa, levels);
  for (std::size_t k = 0; k < a_half_pa.size(); ++k) {
    // Scale-aware: A is compared against the model top, B against unity.
    if (std::abs(a_half_pa[k] - expected.a_half_pa[k]) > 1.0e-9 * top_pressure_pa ||
        std::abs(b_half[k] - expected.b_half[k]) > 1.0e-9) {
      return false;
    }
  }
  return true;
}

void HybridPressureCoefficients::validate(
    const Real minimum_surface_pressure_pa, const Real maximum_surface_pressure_pa,
    const Real minimum_pressure_thickness_pa) const {
  require_positive(minimum_surface_pressure_pa, "minimum surface pressure");
  require_positive(maximum_surface_pressure_pa, "maximum surface pressure");
  require_positive(minimum_pressure_thickness_pa, "minimum pressure thickness");
  if (minimum_surface_pressure_pa > maximum_surface_pressure_pa) {
    throw std::invalid_argument("surface pressure range is inverted");
  }
  if (a_half_pa.size() < 2 || a_half_pa.size() != b_half.size()) {
    throw std::invalid_argument(
        "hybrid coefficient arrays must have equal length >= 2");
  }
  for (std::size_t k = 0; k < a_half_pa.size(); ++k) {
    require_finite(a_half_pa[k], "hybrid A coefficient");
    require_finite(b_half[k], "hybrid B coefficient");
  }
  for (std::size_t k = 1; k < a_half_pa.size(); ++k) {
    for (const Real ps : {minimum_surface_pressure_pa, maximum_surface_pressure_pa}) {
      const Real upper = a_half_pa[k - 1] + b_half[k - 1] * ps;
      const Real lower = a_half_pa[k] + b_half[k] * ps;
      if (!(upper > 0.0) || !(lower - upper >= minimum_pressure_thickness_pa)) {
        throw std::invalid_argument(
            "hybrid interfaces are nonpositive, nonmonotone, or too thin");
      }
    }
  }
}

AtmosphericHybridCoordinate::AtmosphericHybridCoordinate(
    HybridPressureCoefficients coefficients, const Real minimum_surface_pressure_pa,
    const Real maximum_surface_pressure_pa, const Real minimum_pressure_thickness_pa)
    : coefficients_(std::move(coefficients)),
      minimum_surface_pressure_pa_(minimum_surface_pressure_pa),
      maximum_surface_pressure_pa_(maximum_surface_pressure_pa),
      minimum_pressure_thickness_pa_(minimum_pressure_thickness_pa) {
  coefficients_.validate(minimum_surface_pressure_pa_, maximum_surface_pressure_pa_,
                         minimum_pressure_thickness_pa_);
  if (!(coefficients_.a_half_pa.front() > 0.0) || coefficients_.b_half.front() != 0.0 ||
      coefficients_.a_half_pa.back() != 0.0 || coefficients_.b_half.back() != 1.0) {
    throw std::invalid_argument("atmospheric hybrid endpoints are invalid");
  }
  for (std::size_t k = 1; k < coefficients_.b_half.size(); ++k) {
    if (coefficients_.b_half[k] < coefficients_.b_half[k - 1]) {
      throw std::invalid_argument("atmospheric B profile must be nondecreasing");
    }
  }
  if (!(minimum_surface_pressure_pa_ > top_pressure_pa())) {
    throw std::invalid_argument("surface-pressure range must lie above model top");
  }
}

std::size_t AtmosphericHybridCoordinate::levels() const noexcept {
  return coefficients_.a_half_pa.size() - 1;
}

Real AtmosphericHybridCoordinate::top_pressure_pa() const noexcept {
  return coefficients_.a_half_pa.front();
}

const HybridPressureCoefficients& AtmosphericHybridCoordinate::coefficients()
    const noexcept {
  return coefficients_;
}

HybridPressureGeometry AtmosphericHybridCoordinate::geometry(
    const Real surface_pressure_pa, const Real gravity_m_s2,
    const Real gas_constant_j_kg_k, const Real heat_capacity_cp_j_kg_k,
    const Real reference_pressure_pa) const {
  HybridPressureGeometry result;
  geometry(surface_pressure_pa, gravity_m_s2, gas_constant_j_kg_k,
           heat_capacity_cp_j_kg_k, reference_pressure_pa, result);
  return result;
}

void AtmosphericHybridCoordinate::geometry(const Real surface_pressure_pa,
                                           const Real gravity_m_s2,
                                           const Real gas_constant_j_kg_k,
                                           const Real heat_capacity_cp_j_kg_k,
                                           const Real reference_pressure_pa,
                                           HybridPressureGeometry& result) const {
  if (surface_pressure_pa < minimum_surface_pressure_pa_ ||
      surface_pressure_pa > maximum_surface_pressure_pa_) {
    throw std::invalid_argument("surface pressure is outside coordinate range");
  }
  require_positive(gravity_m_s2, "gravity");
  require_positive(gas_constant_j_kg_k, "gas constant");
  require_positive(heat_capacity_cp_j_kg_k, "heat capacity");
  require_positive(reference_pressure_pa, "reference pressure");
  if (!(heat_capacity_cp_j_kg_k > gas_constant_j_kg_k)) {
    throw std::invalid_argument("heat capacity must exceed gas constant");
  }
  const Real kappa = gas_constant_j_kg_k / heat_capacity_cp_j_kg_k;
  result.pressure_half_pa.resize(levels() + 1);
  result.exner_half.resize(levels() + 1);
  result.pressure_full_pa.resize(levels());
  result.exner_full.resize(levels());
  result.delta_pressure_pa.resize(levels());
  result.air_mass_kg_m2.resize(levels());
  for (std::size_t k = 0; k <= levels(); ++k) {
    const Real pressure =
        coefficients_.a_half_pa[k] + coefficients_.b_half[k] * surface_pressure_pa;
    if (!(pressure > 0.0) || !std::isfinite(pressure)) {
      throw std::runtime_error("hybrid interface pressure is invalid");
    }
    result.pressure_half_pa[k] = pressure;
    result.exner_half[k] = std::pow(pressure / reference_pressure_pa, kappa);
  }
  for (std::size_t k = 0; k < levels(); ++k) {
    const Real thickness = result.pressure_half_pa[k + 1] - result.pressure_half_pa[k];
    if (!(thickness >= minimum_pressure_thickness_pa_)) {
      throw std::runtime_error("hybrid layer pressure thickness is invalid");
    }
    result.delta_pressure_pa[k] = thickness;
    result.air_mass_kg_m2[k] = thickness / gravity_m_s2;
    result.exner_full[k] =
        stable_exner_mean(result.pressure_half_pa[k], result.pressure_half_pa[k + 1],
                          kappa, reference_pressure_pa);
    result.pressure_full_pa[k] =
        reference_pressure_pa * std::pow(result.exner_full[k], 1.0 / kappa);
    if (!(result.pressure_full_pa[k] > result.pressure_half_pa[k]) ||
        !(result.pressure_full_pa[k] < result.pressure_half_pa[k + 1])) {
      throw std::runtime_error("hybrid full-level pressure is outside layer");
    }
  }
}

}  // namespace mps
