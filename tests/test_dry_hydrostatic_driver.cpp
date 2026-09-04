#include <algorithm>
#include <cmath>
#include <vector>

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

mps::ExperimentConfig held_suarez_config() {
  auto result = config();
  result.planet = {6371220, 7.29212e-5, 9.80616, 287, 1004, 100000};
  result.run = {0, 20, 10, 0};
  result.dry_hydrostatic.test_case = mps::DryHydrostaticTestCase::kHeldSuarez;
  result.physics.kind = mps::PhysicsKind::kHeldSuarez;
  return result;
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

MPS_TEST_CASE("Held-Suarez forcing is part of every accepted SSP-RK3 step") {
  const auto parameters = held_suarez_config();
  mps::DryHydrostaticDriver driver(parameters);
  auto state = driver.initial_state();
  const auto initial_mass = driver.diagnose(state).air_mass_kg_m2;
  const auto initial_tracer = state.tracer_mass_kg_m2;
  std::vector<mps::DryHydrostaticStepDiagnostics> samples;
  driver.advance(
      state, 10.0,
      [&samples](const mps::DryHydrostaticState&, const mps::DryHydrostaticDerived*,
                 const mps::DryHydrostaticStepDiagnostics& diagnostics) {
        samples.push_back(diagnostics);
      });
  MPS_CHECK(samples.size() >= 2);
  MPS_CHECK(samples.front().physics_rates.thermal_energy_rate_w != 0.0);
  MPS_CHECK(samples.back().thermal_energy_contribution_j != 0.0);
  MPS_CHECK(samples.back().rayleigh_drag_energy_contribution_j <= 0.0);
  MPS_CHECK(state.tracer_mass_kg_m2 == initial_tracer);
  const auto final_mass = driver.diagnose(state).air_mass_kg_m2;
  MPS_CHECK_EQ(final_mass.size(), initial_mass.size());
  mps::Real initial_total = 0.0;
  mps::Real final_total = 0.0;
  for (std::size_t index = 0; index < initial_mass.size(); ++index) {
    const auto cell = index / static_cast<std::size_t>(parameters.vertical.levels);
    initial_total += driver.grid().cells()[cell].area_m2 * initial_mass[index];
    final_total += driver.grid().cells()[cell].area_m2 * final_mass[index];
  }
  MPS_CHECK_NEAR(final_total, initial_total, 2e-15 * initial_total);
}

MPS_TEST_CASE("observer intervals only reduce derived diagnostic samples") {
  struct Trace {
    std::size_t calls = 0;
    std::vector<std::uint64_t> sampled_steps;
    mps::Real thermal_energy_j = 0.0;
    mps::Real drag_energy_j = 0.0;
  } every_step, interval;

  const auto run = [](const std::uint64_t interval_steps, Trace& trace) {
    auto parameters = held_suarez_config();
    parameters.diagnostics.interval_steps = interval_steps;
    const mps::DryHydrostaticDriver driver(parameters);
    auto state = driver.initial_state();
    driver.advance(
        state, parameters.run.end_time_s,
        [&trace](const mps::DryHydrostaticState& sampled,
                 const mps::DryHydrostaticDerived* derived,
                 const mps::DryHydrostaticStepDiagnostics& diagnostics) {
          ++trace.calls;
          trace.thermal_energy_j += diagnostics.thermal_energy_contribution_j;
          trace.drag_energy_j += diagnostics.rayleigh_drag_energy_contribution_j;
          if (derived != nullptr) trace.sampled_steps.push_back(sampled.step);
        });
    return state;
  };

  const auto dense_state = run(1, every_step);
  const auto sparse_state = run(3, interval);
  MPS_CHECK(mps::flatten_dry_hydrostatic_state(dense_state, 2) ==
            mps::flatten_dry_hydrostatic_state(sparse_state, 2));
  MPS_CHECK_EQ(every_step.calls, interval.calls);
  MPS_CHECK_EQ(every_step.thermal_energy_j, interval.thermal_energy_j);
  MPS_CHECK_EQ(every_step.drag_energy_j, interval.drag_energy_j);
  MPS_CHECK(interval.sampled_steps.size() < every_step.sampled_steps.size());
  MPS_CHECK_EQ(interval.sampled_steps.front(), 0U);
  MPS_CHECK_EQ(interval.sampled_steps.back(), sparse_state.step);
  for (std::size_t index = 1; index + 1 < interval.sampled_steps.size(); ++index)
    MPS_CHECK_EQ(interval.sampled_steps[index] % 3, 0U);
}

MPS_TEST_CASE("none physics leaves the dry RHS exactly unchanged") {
  auto unforced_parameters = held_suarez_config();
  unforced_parameters.physics.kind = mps::PhysicsKind::kNone;
  mps::DryHydrostaticDriver unforced(unforced_parameters);
  const auto state = unforced.initial_state();
  const auto before = unforced.rhs(state);

  auto explicit_none = unforced_parameters;
  explicit_none.physics = {.kind = mps::PhysicsKind::kNone};
  mps::DryHydrostaticDriver repeated(explicit_none);
  const auto after = repeated.rhs(state);
  MPS_CHECK(before.surface_pressure_pa_s == after.surface_pressure_pa_s);
  MPS_CHECK(before.tendency.air_mass == after.tendency.air_mass);
  MPS_CHECK(before.tendency.potential_temperature_mass ==
            after.tendency.potential_temperature_mass);
  MPS_CHECK(before.tendency.tracer_mass == after.tendency.tracer_mass);
  MPS_CHECK_EQ(before.physics_diagnostics.thermal_energy_rate_w, 0.0);
  MPS_CHECK_EQ(after.physics_diagnostics.rayleigh_drag_work_w, 0.0);
  for (std::size_t index = 0; index < before.tendency.momentum.size(); ++index) {
    MPS_CHECK_EQ(before.tendency.momentum[index].x, after.tendency.momentum[index].x);
    MPS_CHECK_EQ(before.tendency.momentum[index].y, after.tendency.momentum[index].y);
    MPS_CHECK_EQ(before.tendency.momentum[index].z, after.tendency.momentum[index].z);
  }
}

MPS_TEST_CASE("source-only SSP-RK3 relaxation has third-order time convergence") {
  const auto rate =
      mps::held_suarez_rates(0.0, 100000.0, 100000.0, held_suarez_config().planet)
          .temperature_relaxation_rate_s_1;
  constexpr mps::Real equilibrium = 315.0;
  constexpr mps::Real initial = 264.0;
  constexpr mps::Real end_time = 86400.0;
  const auto integrate = [&](const int steps) {
    const mps::Real dt = end_time / static_cast<mps::Real>(steps);
    mps::Real value = initial;
    const auto rhs = [&](const mps::Real temperature) {
      return -rate * (temperature - equilibrium);
    };
    for (int step = 0; step < steps; ++step) {
      const mps::Real stage1 = value + dt * rhs(value);
      const mps::Real stage2 = 0.75 * value + 0.25 * (stage1 + dt * rhs(stage1));
      value = value / 3.0 + 2.0 / 3.0 * (stage2 + dt * rhs(stage2));
    }
    return value;
  };
  const mps::Real exact =
      equilibrium + (initial - equilibrium) * std::exp(-rate * end_time);
  const mps::Real coarse_error = std::abs(integrate(2) - exact);
  const mps::Real medium_error = std::abs(integrate(4) - exact);
  const mps::Real fine_error = std::abs(integrate(8) - exact);
  MPS_CHECK(coarse_error / medium_error > 7.0);
  MPS_CHECK(medium_error / fine_error > 7.0);
}
int main() { return mps::test::run_all(); }
