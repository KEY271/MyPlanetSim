#include <algorithm>
#include <cmath>

#include "myplanetsim/dynamics/dry_hydrostatic_driver.hpp"
#include "support/test.hpp"

namespace {

mps::ExperimentConfig config() {
  return {.kind = mps::ExperimentKind::kDryHydrostatic,
          .planet = {2, 0, 10, 287, 1004, 100000},
          .run = {0, 1, .1, 0},
          .grid = {2},
          .vertical = {.levels = 2,
                       .a_half_pa = {1000, 500, 0},
                       .b_half = {0, .5, 1},
                       .surface_pressure_pa = 100000,
                       .minimum_surface_pressure_pa = 90000,
                       .maximum_surface_pressure_pa = 110000,
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
}

mps::DryHydrostaticState update(const mps::DryHydrostaticState& base,
                                const mps::DryHydrostaticRhs& rhs, const double scale) {
  auto result = base;
  for (std::size_t c = 0; c < result.surface_pressure_pa.size(); ++c)
    result.surface_pressure_pa[c] += scale * rhs.surface_pressure_pa_s[c];
  for (std::size_t n = 0; n < result.horizontal_momentum_mass_kg_m_s.size(); ++n) {
    result.horizontal_momentum_mass_kg_m_s[n] =
        result.horizontal_momentum_mass_kg_m_s[n] + scale * rhs.tendency.momentum[n];
    result.potential_temperature_mass_k_kg_m2[n] +=
        scale * rhs.tendency.potential_temperature_mass[n];
    result.tracer_mass_kg_m2[n] += scale * rhs.tendency.tracer_mass[n];
  }
  return result;
}

mps::DryHydrostaticState combine(const mps::DryHydrostaticState& first,
                                 const double first_weight,
                                 const mps::DryHydrostaticState& second,
                                 const double second_weight) {
  auto result = first;
  for (std::size_t c = 0; c < result.surface_pressure_pa.size(); ++c)
    result.surface_pressure_pa[c] = first_weight * first.surface_pressure_pa[c] +
                                    second_weight * second.surface_pressure_pa[c];
  for (std::size_t n = 0; n < result.horizontal_momentum_mass_kg_m_s.size(); ++n) {
    result.horizontal_momentum_mass_kg_m_s[n] =
        first_weight * first.horizontal_momentum_mass_kg_m_s[n] +
        second_weight * second.horizontal_momentum_mass_kg_m_s[n];
    result.potential_temperature_mass_k_kg_m2[n] =
        first_weight * first.potential_temperature_mass_k_kg_m2[n] +
        second_weight * second.potential_temperature_mass_k_kg_m2[n];
    result.tracer_mass_kg_m2[n] = first_weight * first.tracer_mass_kg_m2[n] +
                                  second_weight * second.tracer_mass_kg_m2[n];
  }
  return result;
}

void project(const mps::DryHydrostaticDriver& driver, mps::DryHydrostaticState& state,
             const std::size_t levels) {
  for (std::size_t n = 0; n < state.horizontal_momentum_mass_kg_m_s.size(); ++n) {
    const auto cell = n / levels;
    state.horizontal_momentum_mass_kg_m_s[n] = mps::project_tangent(
        state.horizontal_momentum_mass_kg_m_s[n], driver.grid().cells()[cell].center);
  }
}

}  // namespace

MPS_TEST_CASE("uniform isothermal rest remains stationary") {
  const auto c = config();
  mps::DryHydrostaticDriver d(c);
  auto s = d.initial_state();
  auto before = mps::flatten_dry_hydrostatic_state(s, 2);
  d.advance(s, .2);
  auto after = mps::flatten_dry_hydrostatic_state(s, 2);
  MPS_CHECK_EQ(before.size(), after.size());
  for (std::size_t i = 0; i < before.size(); ++i)
    MPS_CHECK_NEAR(before[i], after[i], 1e-10);
  MPS_CHECK(s.step > 0);
}

MPS_TEST_CASE("advance limits the step to the stable CFL and lands on the end time") {
  const auto c = config();
  mps::DryHydrostaticDriver driver(c);
  auto state = driver.initial_state();
  const auto rhs = driver.rhs(state);
  // The Rusanov external-mode estimate is well below the requested maximum step, so the
  // run below is CFL limited rather than limited by run.time_step_s.
  MPS_CHECK(rhs.horizontal_stable_time_step_s < c.run.time_step_s);
  MPS_CHECK(rhs.vertical_stable_time_step_s > 0.0);

  constexpr double end = 0.05;
  driver.advance(state, end);
  MPS_CHECK_EQ(state.time_s, end);
  MPS_CHECK(state.step > 1);
  for (const auto pressure : state.surface_pressure_pa) {
    MPS_CHECK(pressure >= c.vertical.minimum_surface_pressure_pa);
    MPS_CHECK(pressure <= c.vertical.maximum_surface_pressure_pa);
  }
}

MPS_TEST_CASE("advance reevaluates the complete RHS at all SSP-RK3 stages") {
  auto c = config();
  constexpr double dt = 1.0e-5;
  c.run.time_step_s = dt;
  mps::DryHydrostaticDriver driver(c);
  auto initial = driver.initial_state();
  const auto derived = driver.diagnose(initial);
  const auto tangent = mps::normalize(
      mps::project_tangent(mps::Vec3{0, 0, 1}, driver.grid().cells()[0].center));
  initial.horizontal_momentum_mass_kg_m_s[0] =
      0.01 * derived.air_mass_kg_m2[0] * tangent;

  const auto rhs1 = driver.rhs(initial);
  auto stage1 = update(initial, rhs1, dt);
  stage1.time_s = dt;
  project(driver, stage1, c.vertical.levels);

  const auto rhs2 = driver.rhs(stage1);
  auto stage2 = combine(initial, 0.75, update(stage1, rhs2, dt), 0.25);
  stage2.time_s = 0.5 * dt;
  project(driver, stage2, c.vertical.levels);

  const auto rhs3 = driver.rhs(stage2);
  auto expected = combine(initial, 1.0 / 3.0, update(stage2, rhs3, dt), 2.0 / 3.0);
  expected.time_s = dt;
  expected.step = 1;
  project(driver, expected, c.vertical.levels);

  auto actual = initial;
  driver.advance(actual, dt);
  const auto expected_values =
      mps::flatten_dry_hydrostatic_state(expected, c.vertical.levels);
  const auto actual_values =
      mps::flatten_dry_hydrostatic_state(actual, c.vertical.levels);
  MPS_CHECK_EQ(actual.step, expected.step);
  MPS_CHECK_NEAR(actual.time_s, expected.time_s, 0.0);
  MPS_CHECK_EQ(actual_values.size(), expected_values.size());
  for (std::size_t n = 0; n < actual_values.size(); ++n) {
    const double tolerance = 2.0e-13 * std::max(1.0, std::abs(expected_values[n]));
    MPS_CHECK_NEAR(actual_values[n], expected_values[n], tolerance);
  }
}
int main() { return mps::test::run_all(); }
