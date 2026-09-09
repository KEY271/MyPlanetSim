#include "myplanetsim/physics/physics_schedule.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace mps {
namespace {

[[nodiscard]] bool same_time(const Real left, const Real right) {
  return std::abs(left - right) <=
         16.0 * std::numeric_limits<Real>::epsilon() *
             std::max({1.0, std::abs(left), std::abs(right)});
}

void add_deadlines(std::vector<Real>& deadlines, const Real interval,
                   const Real end) {
  for (std::size_t count = 1;; ++count) {
    const Real time = std::min(end, static_cast<Real>(count) * interval);
    deadlines.push_back(time);
    if (same_time(time, end)) break;
  }
}

[[nodiscard]] bool is_due(const Real time, const Real interval, const Real end) {
  if (same_time(time, end)) return true;
  const Real multiple = std::round(time / interval);
  return multiple >= 1.0 && same_time(time, multiple * interval);
}

}  // namespace

std::vector<PhysicsEvent> make_physics_events(
    const PhysicsScheduleParameters& schedule, const Real dynamics_interval_s,
    const bool radiation_active, const bool boundary_layer_active,
    const bool convection_active, const bool moisture_active) {
  if (!(dynamics_interval_s > 0.0) || !std::isfinite(dynamics_interval_s))
    throw std::invalid_argument("physics dynamics interval must be finite and positive");

  std::vector<Real> deadlines;
  if (schedule.kind == PhysicsScheduleKind::kLegacy) {
    deadlines.push_back(dynamics_interval_s);
  } else {
    if (radiation_active)
      add_deadlines(deadlines, schedule.radiation_diagnostic_interval_s,
                    dynamics_interval_s);
    if (boundary_layer_active)
      add_deadlines(deadlines, schedule.boundary_layer_maximum_update_interval_s,
                    dynamics_interval_s);
    if (convection_active)
      add_deadlines(deadlines, schedule.convection_diagnostic_interval_s,
                    dynamics_interval_s);
    if (deadlines.empty()) deadlines.push_back(dynamics_interval_s);
    std::ranges::sort(deadlines);
    const auto unique_end = std::ranges::unique(deadlines, same_time);
    deadlines.erase(unique_end.begin(), unique_end.end());
  }

  std::vector<PhysicsEvent> result;
  result.reserve(deadlines.size());
  Real previous = 0.0;
  Real previous_radiation = 0.0;
  Real previous_boundary_layer = 0.0;
  Real previous_convection = 0.0;
  for (const Real time : deadlines) {
    const bool legacy = schedule.kind == PhysicsScheduleKind::kLegacy;
    const bool radiation = radiation_active &&
                           (legacy || is_due(time, schedule.radiation_diagnostic_interval_s,
                                             dynamics_interval_s));
    const bool boundary_layer =
        boundary_layer_active &&
        (legacy || is_due(time, schedule.boundary_layer_maximum_update_interval_s,
                          dynamics_interval_s));
    const bool convection =
        convection_active &&
        (legacy || is_due(time, schedule.convection_diagnostic_interval_s,
                          dynamics_interval_s));
    result.push_back({.end_time_offset_s = time,
                      .interval_s = time - previous,
                      .radiation_interval_s =
                          radiation ? time - previous_radiation : 0.0,
                      .boundary_layer_interval_s =
                          boundary_layer ? time - previous_boundary_layer : 0.0,
                      .convection_interval_s =
                          convection ? time - previous_convection : 0.0,
                      .radiation = radiation,
                      .boundary_layer = boundary_layer,
                      .convection = convection,
                      .saturation_adjustment =
                          moisture_active && (radiation || boundary_layer || convection)});
    previous = time;
    if (radiation) previous_radiation = time;
    if (boundary_layer) previous_boundary_layer = time;
    if (convection) previous_convection = time;
  }
  return result;
}

}  // namespace mps
