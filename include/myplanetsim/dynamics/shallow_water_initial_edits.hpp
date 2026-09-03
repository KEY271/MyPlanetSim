#pragma once

#include <span>

#include "myplanetsim/control/control_request.hpp"
#include "myplanetsim/dynamics/shallow_water_state.hpp"

namespace mps {

struct InitialConditionEditDiagnostics {
  Real added_volume_m3 = 0.0;
  Real minimum_depth_m = 0.0;
  Real maximum_depth_m = 0.0;
};

void validate_initial_condition_edit(const InitialConditionEditV1& edit);
[[nodiscard]] InitialConditionEditDiagnostics apply_initial_condition_edits(
    const CubedSphereGrid& grid, ShallowWaterState& state,
    std::span<const InitialConditionEditV1> edits, Real depth_floor_m);

}  // namespace mps
