#include "myplanetsim/vertical/vertical_transport.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "myplanetsim/core/validation.hpp"

namespace mps {

VerticalMassFlux diagnose_vertical_mass_flux(
    const std::vector<Real>& horizontal_air_mass_tendency_kg_m2_s,
    const std::vector<Real>& b_half, const Real gravity_m_s2) {
  VerticalMassFlux result;
  diagnose_vertical_mass_flux(horizontal_air_mass_tendency_kg_m2_s, b_half,
                              gravity_m_s2, result);
  return result;
}

void diagnose_vertical_mass_flux(
    const std::span<const Real> horizontal_air_mass_tendency_kg_m2_s,
    const std::span<const Real> b_half, const Real gravity_m_s2,
    VerticalMassFlux& result) {
  const std::size_t nz = horizontal_air_mass_tendency_kg_m2_s.size();
  if (nz == 0 || b_half.size() != nz + 1) {
    throw std::invalid_argument("mass-flux tendency and B shapes differ");
  }
  require_positive(gravity_m_s2, "gravity");
  result.target_air_mass_tendency_kg_m2_s.resize(nz);
  result.interface_flux_kg_m2_s.assign(nz + 1, 0.0);
  Real sum = 0.0;
  Real absolute_sum = 0.0;
  for (const Real tendency : horizontal_air_mass_tendency_kg_m2_s) {
    require_finite(tendency, "horizontal air-mass tendency");
    sum += tendency;
    absolute_sum += std::abs(tendency);
  }
  result.surface_pressure_tendency_pa_s = gravity_m_s2 * sum;
  for (std::size_t k = 0; k < nz; ++k) {
    require_finite(b_half[k], "hybrid B coefficient");
    require_finite(b_half[k + 1], "hybrid B coefficient");
    result.target_air_mass_tendency_kg_m2_s[k] = (b_half[k + 1] - b_half[k]) *
                                                 result.surface_pressure_tendency_pa_s /
                                                 gravity_m_s2;
    result.interface_flux_kg_m2_s[k + 1] = result.interface_flux_kg_m2_s[k] +
                                           horizontal_air_mass_tendency_kg_m2_s[k] -
                                           result.target_air_mass_tendency_kg_m2_s[k];
  }
  result.continuity_residual_pa_s = gravity_m_s2 * result.interface_flux_kg_m2_s.back();
  const Real residual_scale =
      std::max({1.0, std::abs(result.surface_pressure_tendency_pa_s),
                gravity_m_s2 * absolute_sum});
  const Real residual_tolerance =
      64.0 * std::numeric_limits<Real>::epsilon() * residual_scale;
  if (std::abs(result.continuity_residual_pa_s) > residual_tolerance) {
    throw std::runtime_error(
        "vertical mass-flux continuity residual exceeds tolerance");
  }
}

Real vertical_stable_time_step(
    const std::span<const Real> air_mass_kg_m2,
    const std::span<const Real> interface_flux_kg_m2_s,
    const std::span<const Real> horizontal_air_mass_tendency_kg_m2_s, const Real cfl,
    const Real maximum_time_step_s) {
  const std::size_t nz = air_mass_kg_m2.size();
  if (nz == 0 || interface_flux_kg_m2_s.size() != nz + 1 ||
      horizontal_air_mass_tendency_kg_m2_s.size() != nz) {
    throw std::invalid_argument("vertical CFL shapes differ");
  }
  require_positive(cfl, "vertical CFL");
  if (cfl > 1.0) {
    throw std::invalid_argument("vertical CFL must not exceed one");
  }
  require_positive(maximum_time_step_s, "maximum time step");
  Real result = maximum_time_step_s;
  for (std::size_t k = 0; k < nz; ++k) {
    require_positive(air_mass_kg_m2[k], "air mass");
    require_finite(interface_flux_kg_m2_s[k], "vertical interface mass flux");
    require_finite(interface_flux_kg_m2_s[k + 1], "vertical interface mass flux");
    require_finite(horizontal_air_mass_tendency_kg_m2_s[k],
                   "horizontal air-mass tendency");
    const Real outward = std::max(-interface_flux_kg_m2_s[k], 0.0) +
                         std::max(interface_flux_kg_m2_s[k + 1], 0.0) +
                         std::max(-horizontal_air_mass_tendency_kg_m2_s[k], 0.0);
    if (outward > 0.0) {
      result = std::min(result, cfl * air_mass_kg_m2[k] / outward);
    }
  }
  if (!(result > 0.0) || !std::isfinite(result)) {
    throw std::runtime_error("vertical stable time step is invalid");
  }
  return result;
}

Real vertical_maximum_cfl(
    const std::span<const Real> air_mass_kg_m2,
    const std::span<const Real> interface_flux_kg_m2_s,
    const std::span<const Real> horizontal_air_mass_tendency_kg_m2_s,
    const Real time_step_s) {
  const std::size_t nz = air_mass_kg_m2.size();
  if (nz == 0 || interface_flux_kg_m2_s.size() != nz + 1 ||
      horizontal_air_mass_tendency_kg_m2_s.size() != nz) {
    throw std::invalid_argument("vertical CFL shapes differ");
  }
  require_positive(time_step_s, "vertical CFL time step");
  Real result = 0.0;
  for (std::size_t k = 0; k < nz; ++k) {
    require_positive(air_mass_kg_m2[k], "air mass");
    require_finite(interface_flux_kg_m2_s[k], "vertical interface mass flux");
    require_finite(interface_flux_kg_m2_s[k + 1], "vertical interface mass flux");
    require_finite(horizontal_air_mass_tendency_kg_m2_s[k],
                   "horizontal air-mass tendency");
    const Real outward = std::max(-interface_flux_kg_m2_s[k], 0.0) +
                         std::max(interface_flux_kg_m2_s[k + 1], 0.0) +
                         std::max(-horizontal_air_mass_tendency_kg_m2_s[k], 0.0);
    result = std::max(result, time_step_s * outward / air_mass_kg_m2[k]);
  }
  if (!std::isfinite(result)) {
    throw std::runtime_error("vertical CFL is not finite");
  }
  return result;
}

std::vector<Real> vertical_scalar_rhs(
    const std::vector<Real>& scalar_mass, const std::vector<Real>& air_mass_kg_m2,
    const std::vector<Real>& horizontal_scalar_tendency,
    const VerticalMassFlux& mass_flux, const VerticalTransportScheme scheme,
    const VerticalLimiterKind limiter) {
  const std::size_t nz = scalar_mass.size();
  if (nz == 0 || air_mass_kg_m2.size() != nz ||
      horizontal_scalar_tendency.size() != nz ||
      mass_flux.interface_flux_kg_m2_s.size() != nz + 1) {
    throw std::invalid_argument("vertical scalar transport shapes differ");
  }
  std::vector<Real> scalar(nz);
  std::vector<Real> interface_coordinate(nz + 1, 0.0);
  std::vector<Real> center_coordinate(nz);
  std::vector<Real> flux(nz + 1);
  std::vector<Real> rhs(nz);
  for (std::size_t k = 0; k < nz; ++k) {
    interface_coordinate[k + 1] = interface_coordinate[k] + air_mass_kg_m2[k];
    center_coordinate[k] =
        0.5 * (interface_coordinate[k] + interface_coordinate[k + 1]);
  }
  vertical_scalar_rhs(scalar_mass, air_mass_kg_m2, horizontal_scalar_tendency,
                      mass_flux, scheme, limiter, interface_coordinate,
                      center_coordinate, scalar, flux, rhs);
  return rhs;
}

void vertical_scalar_rhs(
    const std::span<const Real> scalar_mass, const std::span<const Real> air_mass_kg_m2,
    const std::span<const Real> horizontal_scalar_tendency,
    const VerticalMassFlux& mass_flux, const VerticalTransportScheme scheme,
    const VerticalLimiterKind limiter, const std::span<const Real> interface_coordinate,
    const std::span<const Real> center_coordinate, const std::span<Real> scalar,
    const std::span<Real> flux, const std::span<Real> rhs) {
  const std::size_t nz = scalar_mass.size();
  if (nz == 0 || air_mass_kg_m2.size() != nz ||
      horizontal_scalar_tendency.size() != nz ||
      mass_flux.interface_flux_kg_m2_s.size() != nz + 1 || scalar.size() != nz ||
      interface_coordinate.size() != nz + 1 || center_coordinate.size() != nz ||
      flux.size() != nz + 1 || rhs.size() != nz) {
    throw std::invalid_argument("vertical scalar transport shapes differ");
  }
  std::fill(flux.begin(), flux.end(), 0.0);
  for (std::size_t k = 0; k < nz; ++k) {
    require_positive(air_mass_kg_m2[k], "air mass");
    require_finite(scalar_mass[k], "scalar mass");
    require_finite(horizontal_scalar_tendency[k], "horizontal scalar tendency");
    scalar[k] = scalar_mass[k] / air_mass_kg_m2[k];
  }
  for (std::size_t interface = 1; interface < nz; ++interface) {
    const Real f = mass_flux.interface_flux_kg_m2_s[interface];
    require_finite(f, "vertical interface mass flux");
    const std::size_t donor = f >= 0.0 ? interface - 1 : interface;
    Real face = scalar[donor];
    if (scheme == VerticalTransportScheme::kLinear && nz > 1) {
      const std::size_t k = donor;
      Real slope = 0.0;
      if (k > 0 && k + 1 < nz) {
        const Real left = (scalar[k] - scalar[k - 1]) /
                          (center_coordinate[k] - center_coordinate[k - 1]);
        const Real right = (scalar[k + 1] - scalar[k]) /
                           (center_coordinate[k + 1] - center_coordinate[k]);
        slope = 0.5 * (left + right);
        if (limiter == VerticalLimiterKind::kMinmod) {
          slope = left * right <= 0.0
                      ? 0.0
                      : std::copysign(std::min(std::abs(left), std::abs(right)), left);
        }
      }
      face += slope * (interface_coordinate[interface] - center_coordinate[k]);
    }
    flux[interface] = f * face;
  }
  for (std::size_t k = 0; k < nz; ++k) {
    rhs[k] = horizontal_scalar_tendency[k] + flux[k] - flux[k + 1];
  }
}

}  // namespace mps
