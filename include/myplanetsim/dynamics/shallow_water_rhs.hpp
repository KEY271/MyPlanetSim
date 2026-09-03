#pragma once

#include "myplanetsim/dynamics/shallow_water_flux.hpp"
#include "myplanetsim/dynamics/shallow_water_state.hpp"

namespace mps {

struct ShallowWaterRhsComponents {
  ShallowWaterTendency flux;
  ShallowWaterTendency pressure;
  ShallowWaterTendency coriolis;
  ShallowWaterTendency diffusion;
  ShallowWaterTendency total;
  Real maximum_wave_speed_m_s;
};

[[nodiscard]] ShallowWaterRhsComponents assemble_first_order_shallow_water_rhs(
    const CubedSphereGrid& grid, const ShallowWaterState& state, Real gravity_m_s2,
    Real rotation_rate_rad_s, Real depth_floor_m);

}  // namespace mps
