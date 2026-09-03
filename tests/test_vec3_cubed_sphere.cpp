#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>

#include "myplanetsim/geometry/cubed_sphere.hpp"
#include "support/test.hpp"

MPS_TEST_CASE("Vec3 provides safe spherical primitives") {
  constexpr mps::Vec3 x{1.0, 0.0, 0.0};
  constexpr mps::Vec3 y{0.0, 1.0, 0.0};
  MPS_CHECK_EQ(mps::dot(x, y), 0.0);
  MPS_CHECK_EQ(mps::cross(x, y).z, 1.0);
  MPS_CHECK_NEAR(mps::safe_angle(x, -x), std::numbers::pi, 2.0e-15);
  MPS_CHECK_THROWS_AS(mps::normalize({}), std::invalid_argument);
  MPS_CHECK_THROWS_AS(
      mps::normalize({std::numeric_limits<double>::quiet_NaN(), 0.0, 0.0}),
      std::invalid_argument);
}

MPS_TEST_CASE("all panel bases are right handed") {
  for (const auto panel : mps::kPanels) {
    const auto basis = mps::panel_basis(panel);
    MPS_CHECK_NEAR(mps::dot(mps::cross(basis.alpha, basis.beta), basis.normal), 1.0,
                   0.0);
  }
}

MPS_TEST_CASE("gnomonic interior mapping round trips") {
  for (const auto panel : mps::kPanels) {
    for (const double alpha : {-0.6, -0.2, 0.0, 0.31, 0.6}) {
      for (const double beta : {-0.55, 0.0, 0.43}) {
        const auto position = mps::map_to_unit_sphere(panel, alpha, beta);
        const auto inverse = mps::inverse_map(position);
        MPS_CHECK(inverse.panel == panel);
        MPS_CHECK_NEAR(inverse.alpha, alpha, 4.0e-16);
        MPS_CHECK_NEAR(inverse.beta, beta, 4.0e-16);
        MPS_CHECK_NEAR(mps::norm(position), 1.0, 4.0e-16);
        const auto tangent = mps::tangent_basis(panel, alpha, beta);
        MPS_CHECK_NEAR(mps::dot(tangent.alpha, position), 0.0, 4.0e-16);
        MPS_CHECK_NEAR(mps::dot(tangent.beta, position), 0.0, 4.0e-16);
        const mps::TangentComponents components{2.0, -3.0};
        const auto round_trip = mps::to_tangent_components(
            mps::from_tangent_components(components, tangent), tangent);
        MPS_CHECK_NEAR(round_trip.alpha, components.alpha, 2.0e-15);
        MPS_CHECK_NEAR(round_trip.beta, components.beta, 2.0e-15);
      }
    }
  }
}

int main() { return mps::test::run_all(); }
