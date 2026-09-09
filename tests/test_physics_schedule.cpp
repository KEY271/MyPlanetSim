#include <stdexcept>

#include "myplanetsim/physics/physics_schedule.hpp"
#include "support/test.hpp"

MPS_TEST_CASE("legacy physics schedule is one ordered event") {
  const auto events = mps::make_physics_events({}, 225.0, true, true, true, true);
  MPS_CHECK_EQ(events.size(), 1U);
  MPS_CHECK_NEAR(events[0].end_time_offset_s, 225.0, 0.0);
  MPS_CHECK_NEAR(events[0].interval_s, 225.0, 0.0);
  MPS_CHECK_NEAR(events[0].boundary_layer_interval_s, 225.0, 0.0);
  MPS_CHECK(events[0].radiation);
  MPS_CHECK(events[0].boundary_layer);
  MPS_CHECK(events[0].convection);
  MPS_CHECK(events[0].saturation_adjustment);
}

MPS_TEST_CASE("process intervals form a deterministic event union") {
  const mps::PhysicsScheduleParameters schedule{
      .kind = mps::PhysicsScheduleKind::kProcessIntervals,
      .boundary_layer_maximum_update_interval_s = 600.0,
      .convection_diagnostic_interval_s = 900.0,
      .convection_update_mode = mps::ConvectionUpdateMode::kIntermittent,
      .radiation_diagnostic_interval_s = 1800.0};
  const auto events =
      mps::make_physics_events(schedule, 1800.0, true, true, true, true);
  MPS_CHECK_EQ(events.size(), 4U);
  MPS_CHECK_NEAR(events[0].end_time_offset_s, 600.0, 0.0);
  MPS_CHECK(events[0].boundary_layer);
  MPS_CHECK(!events[0].convection);
  MPS_CHECK(!events[0].radiation);
  MPS_CHECK_NEAR(events[1].end_time_offset_s, 900.0, 0.0);
  MPS_CHECK_NEAR(events[1].interval_s, 300.0, 0.0);
  MPS_CHECK(events[1].convection);
  MPS_CHECK_NEAR(events[2].end_time_offset_s, 1200.0, 0.0);
  MPS_CHECK(events[2].boundary_layer);
  MPS_CHECK_NEAR(events[2].boundary_layer_interval_s, 600.0, 0.0);
  MPS_CHECK_NEAR(events[3].end_time_offset_s, 1800.0, 0.0);
  MPS_CHECK(events[3].boundary_layer);
  MPS_CHECK(events[3].convection);
  MPS_CHECK(events[3].radiation);
  MPS_CHECK_NEAR(events[3].convection_interval_s, 900.0, 0.0);
  MPS_CHECK_NEAR(events[3].radiation_interval_s, 1800.0, 0.0);
}

MPS_TEST_CASE("process intervals terminate at a short dynamics boundary") {
  const mps::PhysicsScheduleParameters schedule{
      .kind = mps::PhysicsScheduleKind::kProcessIntervals,
      .boundary_layer_maximum_update_interval_s = 600.0,
      .convection_diagnostic_interval_s = 900.0,
      .radiation_diagnostic_interval_s = 1800.0};
  const auto events = mps::make_physics_events(schedule, 450.0, true, true, true, true);
  MPS_CHECK_EQ(events.size(), 1U);
  MPS_CHECK_NEAR(events[0].interval_s, 450.0, 0.0);
  MPS_CHECK(events[0].boundary_layer);
  MPS_CHECK(events[0].convection);
  MPS_CHECK(events[0].radiation);
  MPS_CHECK_THROWS_AS(mps::make_physics_events(schedule, 0.0, true, true, true, true),
                      std::invalid_argument);
}

int main() { return mps::test::run_all(); }
