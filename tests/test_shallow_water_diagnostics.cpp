#include <algorithm>
#include <cmath>
#include <numbers>
#include <vector>

#include "myplanetsim/diagnostics/shallow_water_diagnostics.hpp"
#include "support/test.hpp"

namespace {

[[nodiscard]] mps::ShallowWaterState rest_state(const mps::CubedSphereGrid& grid,
                                                const mps::Real depth) {
  return {.time_s = 0.0,
          .step = 0,
          .depth = std::vector<mps::Real>(grid.cell_count(), depth),
          .momentum = std::vector<mps::Vec3>(grid.cell_count())};
}

[[nodiscard]] mps::ShallowWaterTendency zero_tendency(
    const mps::CubedSphereGrid& grid) {
  return {.depth = std::vector<mps::Real>(grid.cell_count()),
          .momentum = std::vector<mps::Vec3>(grid.cell_count())};
}

}  // namespace

MPS_TEST_CASE("shallow-water diagnostics recover constant-rest invariants") {
  constexpr mps::Real radius = 2.0;
  constexpr mps::Real depth = 3.0;
  constexpr mps::Real gravity = 4.0;
  constexpr mps::Real omega = 0.1;
  const mps::CubedSphereGrid grid(32, radius);
  const auto state = rest_state(grid, depth);
  const auto diagnostics =
      mps::diagnostics::diagnose_shallow_water(grid, state, gravity, {0.0, 0.0, omega});
  const mps::Real area = grid.total_area_m2();
  MPS_CHECK_NEAR(diagnostics.mass, depth * area, 2.0e-14 * depth * area);
  MPS_CHECK_NEAR(diagnostics.energy, 0.5 * gravity * depth * depth * area,
                 2.0e-14 * diagnostics.energy);
  MPS_CHECK_NEAR(diagnostics.minimum_depth, depth, 0.0);
  MPS_CHECK_NEAR(diagnostics.maximum_depth, depth, 0.0);
  MPS_CHECK_NEAR(diagnostics.maximum_radial_momentum, 0.0, 0.0);
  const mps::Real exact_aam = (2.0 / 3.0) * omega * radius * radius * depth * area;
  MPS_CHECK_NEAR(diagnostics.axial_angular_momentum, exact_aam, 4.0e-4 * exact_aam);
}

MPS_TEST_CASE("shallow-water invariant rates have physical unit scaling") {
  const mps::CubedSphereGrid grid(8, 2.0);
  constexpr mps::Real depth = 3.0;
  constexpr mps::Real gravity = 4.0;
  const auto state = rest_state(grid, depth);
  auto tendency = zero_tendency(grid);
  std::ranges::fill(tendency.depth, 0.25);
  const auto rates = mps::diagnostics::shallow_water_invariant_rates(
      grid, state, tendency, gravity, {});
  MPS_CHECK_NEAR(rates.mass, 0.25 * grid.total_area_m2(),
                 2.0e-14 * grid.total_area_m2());
  MPS_CHECK_NEAR(rates.energy, gravity * depth * 0.25 * grid.total_area_m2(),
                 2.0e-14 * gravity * depth * grid.total_area_m2());
}

MPS_TEST_CASE("shallow-water budget components sum to total") {
  const mps::CubedSphereGrid grid(4, 1.0);
  const auto state = rest_state(grid, 2.0);
  auto flux = zero_tendency(grid);
  auto coriolis = zero_tendency(grid);
  auto pressure = zero_tendency(grid);
  auto diffusion = zero_tendency(grid);
  std::ranges::fill(flux.depth, 0.1);
  std::ranges::fill(diffusion.depth, -0.1);
  const auto budget = mps::diagnostics::make_shallow_water_budget(
      grid, state, flux, coriolis, pressure, diffusion, 3.0, {0.0, 0.0, 0.2});
  MPS_CHECK_NEAR(budget.total.mass, 0.0, 1.0e-14);
  MPS_CHECK_NEAR(budget.residual.mass, 0.0, 1.0e-14);
  MPS_CHECK_NEAR(budget.residual.energy, 0.0, 1.0e-13);
}

int main() { return mps::test::run_all(); }
