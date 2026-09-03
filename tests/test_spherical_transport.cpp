#include <algorithm>
#include <bit>
#include <cmath>
#include <iostream>
#include <numbers>
#include <vector>

#include "myplanetsim/numerics/spherical_operators.hpp"
#include "myplanetsim/transport/spherical_transport.hpp"
#include "support/test.hpp"

namespace {

[[nodiscard]] mps::ExperimentConfig transport_config(const mps::Index n) {
  return {.kind = mps::ExperimentKind::kSphereTransport,
          .planet = {.radius_m = 1.0,
                     .rotation_rate_rad_s = 1.0,
                     .gravity_m_s2 = 1.0,
                     .gas_constant_j_kg_k = 1.0,
                     .heat_capacity_cp_j_kg_k = 2.0,
                     .reference_pressure_pa = 1.0},
          .run = {.start_time_s = 0.0,
                  .end_time_s = 2.0 * std::numbers::pi,
                  .time_step_s = 1.0,
                  .random_seed = 1},
          .ode = {},
          .grid = {.cells_per_panel = n},
          .transport = {.test_case = mps::TransportTestCase::kSolidBody,
                        .initial_condition = mps::InitialConditionKind::kGaussianHill,
                        .scheme = mps::TransportScheme::kUpwind,
                        .limiter = mps::LimiterKind::kNone,
                        .cfl = 0.35,
                        .rotation_axis_x = 0.0,
                        .rotation_axis_y = 0.0,
                        .rotation_axis_z = 1.0,
                        .angular_speed_rad_s = 1.0},
          .output_directory = "unused"};
}

[[nodiscard]] double convergence_order(const double coarse, const double fine) {
  return std::log(coarse / fine) / std::log(2.0);
}

}  // namespace

MPS_TEST_CASE("solid body edge flux is discretely nondivergent") {
  const mps::CubedSphereGrid grid(12, 2.0);
  const auto flux = mps::prescribed_edge_fluxes(
      grid, mps::TransportTestCase::kSolidBody, 0.0, 1.0, {1.0, 2.0, 3.0}, 0.75);
  const auto divergence = mps::finite_volume_divergence(grid, flux);
  for (const double value : divergence) {
    MPS_CHECK_NEAR(value, 0.0, 2.0e-14);
  }
}

MPS_TEST_CASE("transport preserves constants mass and limiter bounds") {
  auto config = transport_config(10);
  config.run.end_time_s = 0.7;
  config.transport.initial_condition = mps::InitialConditionKind::kConstant;
  config.transport.scheme = mps::TransportScheme::kLinear;
  config.transport.limiter = mps::LimiterKind::kBarthJespersen;
  const auto constant = mps::run_spherical_transport(config);
  for (const double value : constant.state.tracer) {
    MPS_CHECK_NEAR(value, 1.0, 3.0e-15);
  }
  MPS_CHECK(std::abs(constant.diagnostics.relative_mass_drift) < 5.0e-13);

  config.transport.initial_condition = mps::InitialConditionKind::kSlottedCylinder;
  const auto bounded = mps::run_spherical_transport(config);
  MPS_CHECK(bounded.diagnostics.minimum >= -1.0e-12);
  MPS_CHECK(bounded.diagnostics.maximum <= 1.0 + 1.0e-12);
  MPS_CHECK(std::abs(bounded.diagnostics.relative_mass_drift) < 5.0e-13);
}

MPS_TEST_CASE("transport restart is bitwise identical") {
  auto config = transport_config(8);
  config.run.end_time_s = 1.3;
  config.transport.scheme = mps::TransportScheme::kLinear;
  config.transport.limiter = mps::LimiterKind::kBarthJespersen;
  const auto continuous = mps::run_spherical_transport(config);
  const auto stopped = mps::run_spherical_transport(config, std::nullopt, 5);
  MPS_CHECK(!stopped.reached_end_time);
  const auto restarted = mps::run_spherical_transport(config, stopped.state);
  MPS_CHECK_EQ(restarted.state.step, continuous.state.step);
  MPS_CHECK_EQ(restarted.state.time_s, continuous.state.time_s);
  MPS_CHECK(restarted.state.tracer == continuous.state.tracer);
}

MPS_TEST_CASE("solid body transport remains conservative for rotated axes") {
  auto config = transport_config(8);
  config.run.end_time_s = 0.8;
  config.transport.scheme = mps::TransportScheme::kLinear;
  config.transport.limiter = mps::LimiterKind::kNone;
  for (const mps::Vec3 axis :
       {mps::Vec3{1.0, 0.0, 0.0}, mps::Vec3{0.0, 0.0, 1.0}, mps::Vec3{1.0, 2.0, 3.0}}) {
    config.transport.rotation_axis_x = axis.x;
    config.transport.rotation_axis_y = axis.y;
    config.transport.rotation_axis_z = axis.z;
    const auto result = mps::run_spherical_transport(config);
    MPS_CHECK(std::isfinite(result.diagnostics.l2_error));
    MPS_CHECK(std::abs(result.diagnostics.relative_mass_drift) < 5.0e-13);
    MPS_CHECK(result.diagnostics.seam_rms_error < 3.0 * result.diagnostics.l2_error);
  }
}

MPS_TEST_CASE("smooth solid body transport reaches planned convergence orders") {
  std::vector<double> first_order_errors;
  std::vector<double> linear_errors;
  for (const mps::Index n : {12, 24, 48}) {
    auto config = transport_config(n);
    auto first_order = mps::run_spherical_transport(config);
    first_order_errors.push_back(first_order.diagnostics.l1_error);
    MPS_CHECK(std::abs(first_order.diagnostics.relative_mass_drift) < 5.0e-13);
  }
  for (const mps::Index n : {6, 12, 24}) {
    auto config = transport_config(n);
    config.transport.scheme = mps::TransportScheme::kLinear;
    config.transport.limiter = mps::LimiterKind::kNone;
    auto linear = mps::run_spherical_transport(config);
    linear_errors.push_back(linear.diagnostics.l2_error);
    MPS_CHECK(std::abs(linear.diagnostics.relative_mass_drift) < 5.0e-13);
  }
  const double first_order =
      convergence_order(first_order_errors[1], first_order_errors[2]);
  const double linear = convergence_order(linear_errors[1], linear_errors[2]);
  std::cout << "first-order errors = " << first_order_errors[0] << ", "
            << first_order_errors[1] << ", " << first_order_errors[2]
            << "; linear errors = " << linear_errors[0] << ", " << linear_errors[1]
            << ", " << linear_errors[2] << "; observed orders = " << first_order << ", "
            << linear << '\n';
  MPS_CHECK(first_order >= 0.8);
  MPS_CHECK(linear >= 1.7);
}

int main() { return mps::test::run_all(); }
