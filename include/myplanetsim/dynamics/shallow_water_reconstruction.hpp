#pragma once

#include <cstdint>
#include <vector>

#include "myplanetsim/config/experiment_config.hpp"
#include "myplanetsim/dynamics/shallow_water_flux.hpp"
#include "myplanetsim/dynamics/shallow_water_state.hpp"

namespace mps {

struct ShallowWaterFaceStates {
  ShallowWaterPrimitive left;
  ShallowWaterPrimitive right;
};

struct ShallowWaterReconstruction {
  std::vector<ShallowWaterFaceStates> edges;
  std::uint64_t limiter_activations;
};

[[nodiscard]] ShallowWaterReconstruction reconstruct_shallow_water_face_states(
    const CubedSphereGrid& grid, const ShallowWaterState& state,
    ReconstructionKind reconstruction, LimiterKind limiter, Real depth_floor_m);

}  // namespace mps
