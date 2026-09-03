#pragma once

#include "myplanetsim/config/experiment_config.hpp"
#include "myplanetsim/dynamics/shallow_water_state.hpp"

namespace mps {

[[nodiscard]] Real stable_diffusion_time_step(const CubedSphereGrid& grid,
                                              DiffusionKind kind,
                                              Real diffusion_coefficient,
                                              Real maximum_time_step_s);
[[nodiscard]] ShallowWaterTendency shallow_water_diffusion_tendency(
    const CubedSphereGrid& grid, const ShallowWaterState& state, DiffusionKind kind,
    Real diffusion_coefficient);

}  // namespace mps
