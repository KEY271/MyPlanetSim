#pragma once

#include <vector>

#include "myplanetsim/config/experiment_config.hpp"

namespace mps {

struct PhysicsEvent {
  Real end_time_offset_s = 0.0;
  Real interval_s = 0.0;
  Real radiation_interval_s = 0.0;
  Real boundary_layer_interval_s = 0.0;
  Real convection_interval_s = 0.0;
  bool radiation = false;
  bool boundary_layer = false;
  bool convection = false;
  bool saturation_adjustment = false;
};

// Construct the event union for one accepted dynamics interval. In process-interval
// mode every active process also runs at the dynamics-step end, so no increment or
// ledger entry is deferred across a retry/restart boundary.
[[nodiscard]] std::vector<PhysicsEvent> make_physics_events(
    const PhysicsScheduleParameters& schedule, Real dynamics_interval_s,
    bool radiation_active, bool boundary_layer_active, bool convection_active,
    bool moisture_active);

}  // namespace mps
