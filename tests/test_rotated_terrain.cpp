#include <algorithm>
#include <array>
#include <cmath>
#include <utility>

#include "myplanetsim/dynamics/surface_orography.hpp"
#include "support/test.hpp"

namespace {

mps::Vec3 rotate(const mps::Vec3 value, const mps::Vec3 axis, const double angle) {
  const auto unit_axis = mps::normalize(axis);
  return std::cos(angle) * value + std::sin(angle) * mps::cross(unit_axis, value) +
         (1.0 - std::cos(angle)) * mps::dot(unit_axis, value) * unit_axis;
}

struct Norms {
  double mean;
  double l2;
  double maximum;
};

Norms terrain_norms(const mps::CubedSphereGrid& grid, const mps::Vec3 axis,
                    const double angle) {
  double integral = 0;
  double square_integral = 0;
  double maximum = 0;
  for (const auto& cell : grid.cells()) {
    const double height =
        mps::dcmip_2_0_0_surface_height_m(rotate(cell.center, axis, -angle));
    integral += cell.area_m2 * height;
    square_integral += cell.area_m2 * height * height;
    maximum = std::max(maximum, height);
  }
  return {.mean = integral / grid.total_area_m2(),
          .l2 = std::sqrt(square_integral / grid.total_area_m2()),
          .maximum = maximum};
}

}  // namespace

MPS_TEST_CASE("rotated smooth terrain has grid-independent integral norms") {
  const mps::CubedSphereGrid grid(24, 6371220);
  const auto baseline = terrain_norms(grid, {0, 0, 1}, 0);
  constexpr std::array rotations{std::pair{mps::Vec3{0, 0, 1}, 0.37},
                                 std::pair{mps::Vec3{1, 1, 0}, 0.61},
                                 std::pair{mps::Vec3{1, -1, 1}, 0.79}};
  for (const auto& [axis, angle] : rotations) {
    const auto rotated = terrain_norms(grid, axis, angle);
    MPS_CHECK_NEAR(rotated.mean, baseline.mean, 0.015 * baseline.mean);
    MPS_CHECK_NEAR(rotated.l2, baseline.l2, 0.015 * baseline.l2);
    MPS_CHECK_NEAR(rotated.maximum, baseline.maximum, 0.03 * baseline.maximum);
  }
}

int main() { return mps::test::run_all(); }
