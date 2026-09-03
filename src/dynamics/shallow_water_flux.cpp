#include "myplanetsim/dynamics/shallow_water_flux.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "myplanetsim/core/validation.hpp"

namespace mps {
namespace {

void validate_primitive(const ShallowWaterPrimitive& state,
                        const EdgeTangentBasis& basis) {
  require_positive(state.depth_m, "Rusanov edge depth");
  if (!is_finite(state.velocity_m_s) || !is_finite(basis.normal) ||
      !is_finite(basis.tangent)) {
    throw std::invalid_argument("Rusanov edge state is non-finite");
  }
  constexpr Real tolerance = 1.0e-12;
  if (std::abs(norm(basis.normal) - 1.0) > tolerance ||
      std::abs(norm(basis.tangent) - 1.0) > tolerance ||
      std::abs(dot(basis.normal, basis.tangent)) > tolerance) {
    throw std::invalid_argument("Rusanov edge basis is not orthonormal");
  }
}

}  // namespace

ShallowWaterEdgeFlux rusanov_shallow_water_flux(const ShallowWaterPrimitive& left,
                                                const ShallowWaterPrimitive& right,
                                                const EdgeTangentBasis& basis,
                                                const Real gravity_m_s2) {
  require_positive(gravity_m_s2, "Rusanov gravity");
  validate_primitive(left, basis);
  validate_primitive(right, basis);
  const Vec3 left_velocity = dot(left.velocity_m_s, basis.normal) * basis.normal +
                             dot(left.velocity_m_s, basis.tangent) * basis.tangent;
  const Vec3 right_velocity = dot(right.velocity_m_s, basis.normal) * basis.normal +
                              dot(right.velocity_m_s, basis.tangent) * basis.tangent;
  const Real left_normal_velocity = dot(left_velocity, basis.normal);
  const Real right_normal_velocity = dot(right_velocity, basis.normal);
  const Real wave_speed = std::max(
      std::abs(left_normal_velocity) + std::sqrt(gravity_m_s2 * left.depth_m),
      std::abs(right_normal_velocity) + std::sqrt(gravity_m_s2 * right.depth_m));
  const Real left_mass_flux = left.depth_m * left_normal_velocity;
  const Real right_mass_flux = right.depth_m * right_normal_velocity;
  const Vec3 left_momentum = left.depth_m * left_velocity;
  const Vec3 right_momentum = right.depth_m * right_velocity;
  const Vec3 left_momentum_flux =
      left_mass_flux * left_velocity +
      0.5 * gravity_m_s2 * left.depth_m * left.depth_m * basis.normal;
  const Vec3 right_momentum_flux =
      right_mass_flux * right_velocity +
      0.5 * gravity_m_s2 * right.depth_m * right.depth_m * basis.normal;
  return {
      .mass_m2_s = 0.5 * (left_mass_flux + right_mass_flux) -
                   0.5 * wave_speed * (right.depth_m - left.depth_m),
      .momentum_m3_s2 = 0.5 * (left_momentum_flux + right_momentum_flux) -
                        0.5 * wave_speed * (right_momentum - left_momentum),
      .maximum_wave_speed_m_s = wave_speed,
  };
}

}  // namespace mps
