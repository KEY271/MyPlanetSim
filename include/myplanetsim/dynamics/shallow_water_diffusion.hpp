#pragma once

#include "myplanetsim/config/experiment_config.hpp"
#include "myplanetsim/dynamics/shallow_water_state.hpp"

namespace mps {

[[nodiscard]] ShallowWaterTendency shallow_water_diffusion_tendency(
    const CubedSphereGrid& grid, const ShallowWaterState& state, DiffusionKind kind,
    Real diffusion_coefficient);

}  // namespace mps
