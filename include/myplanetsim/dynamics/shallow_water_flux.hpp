#pragma once

#include "myplanetsim/numerics/spherical_operators.hpp"

namespace mps {

struct ShallowWaterPrimitive {
  Real depth_m;
  Vec3 velocity_m_s;
};

struct ShallowWaterEdgeFlux {
  Real mass_m2_s;
  Vec3 momentum_m3_s2;
  Real maximum_wave_speed_m_s;
};

[[nodiscard]] ShallowWaterEdgeFlux rusanov_shallow_water_flux(
    const ShallowWaterPrimitive& left, const ShallowWaterPrimitive& right,
    const EdgeTangentBasis& basis, Real gravity_m_s2);

}  // namespace mps
