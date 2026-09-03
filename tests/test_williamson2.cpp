#include <algorithm>
#include <cmath>

#include "myplanetsim/dynamics/shallow_water_benchmarks.hpp"
#include "support/test.hpp"

namespace {

constexpr mps::Real kRadius = 6.37122e6;
constexpr mps::Real kGravity = 9.80616;
constexpr mps::Real kRotation = 7.292e-5;
constexpr mps::Real kDepth = 2.94e4 / kGravity;
constexpr mps::Real kVelocity =
    2.0 * 3.14159265358979323846 * kRadius / (12.0 * 86400.0);

}  // namespace

MPS_TEST_CASE("Williamson 2 has analytic equatorial speed and polar depth") {
  const mps::CubedSphereGrid grid(12, kRadius);
  const auto state = mps::make_williamson2_state(grid, 0.0, kGravity, kRotation, kDepth,
                                                 kVelocity, {0.0, 0.0, 1.0});
  mps::Real maximum_speed = 0.0;
  mps::Real minimum_depth = kDepth;
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
    maximum_speed = std::max(maximum_speed, mps::norm(state.velocity(cell)));
    minimum_depth = std::min(minimum_depth, state.depth[cell]);
    MPS_CHECK_NEAR(mps::dot(state.momentum[cell], grid.cells()[cell].center), 0.0,
                   1.0e-10);
  }
  const mps::Real coefficient =
      (kRadius * kRotation * kVelocity + 0.5 * kVelocity * kVelocity) / kGravity;
  MPS_CHECK(maximum_speed > 0.99 * kVelocity);
  MPS_CHECK(minimum_depth < kDepth - 0.98 * coefficient);
}

MPS_TEST_CASE("Williamson 2 rotates covariantly with its flow axis") {
  const mps::CubedSphereGrid grid(8, kRadius);
  const auto z_axis = mps::make_williamson2_state(grid, 0.0, kGravity, kRotation,
                                                  kDepth, kVelocity, {0.0, 0.0, 1.0});
  const auto x_axis = mps::make_williamson2_state(grid, 0.0, kGravity, kRotation,
                                                  kDepth, kVelocity, {1.0, 0.0, 0.0});
  mps::Real z_mass = 0.0;
  mps::Real x_mass = 0.0;
  mps::Real z_kinetic = 0.0;
  mps::Real x_kinetic = 0.0;
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
    const mps::Real area = grid.cells()[cell].area_m2;
    z_mass += area * z_axis.depth[cell];
    x_mass += area * x_axis.depth[cell];
    z_kinetic += area * z_axis.depth[cell] * mps::norm_squared(z_axis.velocity(cell));
    x_kinetic += area * x_axis.depth[cell] * mps::norm_squared(x_axis.velocity(cell));
  }
  MPS_CHECK_NEAR(x_mass, z_mass, 2.0e-15 * z_mass);
  MPS_CHECK_NEAR(x_kinetic, z_kinetic, 2.0e-15 * z_kinetic);
}

MPS_TEST_CASE("Williamson 2 rejects invalid physical constants") {
  const mps::CubedSphereGrid grid(2, 1.0);
  MPS_CHECK_THROWS_AS(
      mps::make_williamson2_state(grid, 0.0, 1.0, 1.0, 1.0, 1.0, {0.0, 0.0, 0.0}),
      std::invalid_argument);
  MPS_CHECK_THROWS_AS(
      mps::make_williamson2_state(grid, 0.0, 1.0, 1.0, 1.0, 2.0, {0.0, 0.0, 1.0}),
      std::invalid_argument);
}

int main() { return mps::test::run_all(); }
