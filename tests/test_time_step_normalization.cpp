#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include "myplanetsim/dynamics/dry_hydrostatic_driver.hpp"
#include "myplanetsim/dynamics/shallow_water_driver.hpp"
#include "myplanetsim/numerics/spherical_operators.hpp"
#include "support/test.hpp"

// Gates for the per-cell Courant condition of ADR 0011,
//
//   dt <= cfl * A_cell / sum_f( lambda_f * L_f )
//
// The dry core used to divide by a single face instead of summing over the four faces,
// which is about four times weaker on a quadrilateral cell: the Phase 7 production
// candidate reported CFL 0.45 while running at an effective cell CFL near 1.02.

namespace {

constexpr mps::Real kCfl = 0.5;

mps::ExperimentConfig config() {
  mps::ExperimentConfig result{
      .kind = mps::ExperimentKind::kDryHydrostatic,
      .planet = {6371220, 7.29212e-5, 9.80616, 287, 1004, 100000},
      .run = {0, 3600, 3600, 0},
      .grid = {4},
      .vertical = {.levels = 4,
                   .a_half_pa = {1000, 750, 500, 250, 0},
                   .b_half = {0, .25, .5, .75, 1},
                   .surface_pressure_pa = 100000,
                   .minimum_surface_pressure_pa = 50000,
                   .maximum_surface_pressure_pa = 150000,
                   .minimum_pressure_thickness_pa = 100,
                   .initial_temperature_k = 280,
                   .initial_potential_temperature_k = 300,
                   .temperature_floor_k = 100,
                   .transport_scheme = mps::VerticalTransportScheme::kDonorCell,
                   .limiter = mps::VerticalLimiterKind::kNone,
                   .cfl = .5},
      .dry_hydrostatic = {},
      .diagnostics = {1},
      .output_directory = "x"};
  result.dry_hydrostatic.test_case = mps::DryHydrostaticTestCase::kHeldSuarez;
  result.dry_hydrostatic.cfl = kCfl;
  return result;
}

// The Courant number the reported step actually achieves, computed independently of the
// driver: the largest per-cell, per-level sum of face wave speed times face length,
// divided by cell area, times the step.
[[nodiscard]] mps::Real achieved_cell_courant(const mps::DryHydrostaticDriver& driver,
                                              const mps::DryHydrostaticState& state,
                                              const mps::Real time_step_s) {
  const auto& grid = driver.grid();
  const auto derived = driver.diagnose(state);
  const mps::Real gamma = 1004.0 / (1004.0 - 287.0);
  std::vector<mps::Real> denominator(derived.cells * derived.levels, 0.0);
  for (const auto& edge : grid.edges()) {
    const auto basis = mps::edge_tangent_basis(edge);
    for (const auto cell :
         {grid.cell_index(edge.left_cell), grid.cell_index(edge.right_cell)}) {
      for (std::size_t k = 0; k < derived.levels; ++k) {
        const auto n = mps::dry_hydrostatic_offset(cell, k, derived.levels);
        const mps::Real sound = std::sqrt(gamma * 287.0 * derived.temperature_k[n]);
        const mps::Real normal_speed =
            std::abs(mps::dot(derived.velocity_m_s[n], basis.normal));
        denominator[n] += edge.length_m * (normal_speed + sound);
      }
    }
  }
  mps::Real courant = 0.0;
  for (std::size_t cell = 0; cell < derived.cells; ++cell) {
    for (std::size_t k = 0; k < derived.levels; ++k) {
      const auto n = mps::dry_hydrostatic_offset(cell, k, derived.levels);
      courant =
          std::max(courant, time_step_s * denominator[n] / grid.cells()[cell].area_m2);
    }
  }
  return courant;
}

}  // namespace

MPS_TEST_CASE("the dry core reports a per-cell, not a per-edge, stable step") {
  const mps::DryHydrostaticDriver driver(config());
  const auto state = driver.initial_state();
  const auto rhs = driver.rhs(state);
  const mps::Real step = rhs.horizontal_stable_time_step_s;
  MPS_CHECK(step > 0.0 && std::isfinite(step));

  // The reported step must land on the requested CFL, not roughly four times past it.
  // The per-edge form gave about 4 * kCfl here.
  MPS_CHECK_NEAR(achieved_cell_courant(driver, state, step), kCfl, 0.02 * kCfl);
}

MPS_TEST_CASE("the shallow-water solver already uses the same per-cell definition") {
  const mps::CubedSphereGrid grid(4, 6371220.0);
  mps::ShallowWaterState state;
  state.depth.assign(grid.cell_count(), 5000.0);
  state.momentum.assign(grid.cell_count(), {});
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
    state.momentum[cell] =
        state.depth[cell] * 20.0 *
        mps::cross(mps::Vec3{0.0, 0.0, 1.0}, grid.cells()[cell].center);
  }
  mps::project_shallow_water_momentum(grid, state);
  const mps::Real step = mps::stable_shallow_water_time_step(
      grid, state, 9.80616, kCfl, std::numeric_limits<mps::Real>::max());
  // stable_shallow_water_time_step and shallow_water_cfl_number are inverses of one
  // per-cell condition; the dry core now shares that condition.
  MPS_CHECK_NEAR(mps::shallow_water_cfl_number(grid, state, 9.80616, step), kCfl,
                 1.0e-12 * kCfl);
}

int main() { return mps::test::run_all(); }
