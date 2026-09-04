#include <algorithm>
#include <cmath>
#include <sstream>
#include <vector>

#include "myplanetsim/diagnostics/climate_statistics.hpp"
#include "myplanetsim/diagnostics/dry_hydrostatic_diagnostics.hpp"
#include "myplanetsim/dynamics/dry_hydrostatic_driver.hpp"
#include "myplanetsim/physics/held_suarez.hpp"
#include "support/test.hpp"

namespace {

[[nodiscard]] mps::ExperimentConfig config() {
  return {.kind = mps::ExperimentKind::kDryHydrostatic,
          .planet = {6371220, 7.29212e-5, 9.80616, 287, 1004, 100000},
          .run = {0, 30, 5, 11},
          .grid = {2},
          .vertical = {.levels = 2,
                       .a_half_pa = {1000, 500, 0},
                       .b_half = {0, .5, 1},
                       .surface_pressure_pa = 100000,
                       .minimum_surface_pressure_pa = 80000,
                       .maximum_surface_pressure_pa = 120000,
                       .minimum_pressure_thickness_pa = 100,
                       .initial_temperature_k = 264,
                       .initial_potential_temperature_k = 300,
                       .temperature_floor_k = 100,
                       .transport_scheme = mps::VerticalTransportScheme::kLinear,
                       .limiter = mps::VerticalLimiterKind::kMinmod,
                       .cfl = .5},
          .dry_hydrostatic = {.test_case = mps::DryHydrostaticTestCase::kHeldSuarez},
          .physics = {.kind = mps::PhysicsKind::kHeldSuarez},
          .diagnostics = {1},
          .output_directory = "x"};
}

[[nodiscard]] mps::Real temperature_departure(const mps::DryHydrostaticDriver& driver,
                                              const mps::DryHydrostaticState& state,
                                              const mps::DryHydrostaticDerived& derived,
                                              const mps::PlanetParameters& planet) {
  mps::Real sum = 0.0;
  mps::Real area_sum = 0.0;
  for (std::size_t cell = 0; cell < derived.cells; ++cell) {
    const auto area = driver.grid().cells()[cell].area_m2;
    const auto latitude = std::asin(driver.grid().cells()[cell].center.z);
    for (std::size_t level = 0; level < derived.levels; ++level) {
      const auto offset = mps::dry_hydrostatic_offset(cell, level, derived.levels);
      const auto equilibrium =
          mps::held_suarez_rates(latitude, derived.pressure_pa[offset],
                                 state.surface_pressure_pa[cell], planet)
              .equilibrium_temperature_k;
      sum += area * std::abs(derived.temperature_k[offset] - equilibrium);
      area_sum += area;
    }
  }
  return sum / area_sum;
}

}  // namespace

MPS_TEST_CASE(
    "forced dry core remains finite and moves temperature toward equilibrium") {
  const auto parameters = config();
  const mps::DryHydrostaticDriver driver(parameters);
  auto state = driver.initial_state();
  const auto initial_derived = driver.diagnose(state);
  const auto initial_health = mps::diagnose_dry_hydrostatic_budgets(
      driver.grid(), state, initial_derived, parameters.planet);
  const auto initial_departure =
      temperature_departure(driver, state, initial_derived, parameters.planet);
  mps::Real thermal_contribution = 0.0;
  mps::Real drag_contribution = 0.0;
  driver.advance(state, parameters.run.end_time_s,
                 [&](const mps::DryHydrostaticState&, const mps::DryHydrostaticDerived&,
                     const mps::DryHydrostaticStepDiagnostics& diagnostics) {
                   thermal_contribution += diagnostics.thermal_energy_contribution_j;
                   drag_contribution += diagnostics.rayleigh_drag_energy_contribution_j;
                 });
  const auto final_derived = driver.diagnose(state);
  const auto final_health = mps::diagnose_dry_hydrostatic_budgets(
      driver.grid(), state, final_derived, parameters.planet);
  MPS_CHECK(temperature_departure(driver, state, final_derived, parameters.planet) <
            initial_departure);
  MPS_CHECK(thermal_contribution != 0.0);
  MPS_CHECK(drag_contribution <= 0.0);
  MPS_CHECK_NEAR(final_health.dry_mass_kg, initial_health.dry_mass_kg,
                 2e-14 * initial_health.dry_mass_kg);
  MPS_CHECK_EQ(final_health.tracer_mass_kg, initial_health.tracer_mass_kg);
  MPS_CHECK(final_health.minimum_temperature_k >=
            parameters.vertical.temperature_floor_k);
  for (const auto pressure : state.surface_pressure_pa) {
    MPS_CHECK(pressure >= parameters.vertical.minimum_surface_pressure_pa);
    MPS_CHECK(pressure <= parameters.vertical.maximum_surface_pressure_pa);
  }
  for (const auto value : mps::flatten_dry_hydrostatic_state(state, 2))
    MPS_CHECK(std::isfinite(value));
}

MPS_TEST_CASE(
    "restart before the statistics window reproduces state budget and moments") {
  const auto parameters = config();
  const mps::DryHydrostaticDriver driver(parameters);
  constexpr mps::Real restart_time = 10.0;
  constexpr mps::Real statistics_start = 20.0;

  auto uninterrupted = driver.initial_state();
  mps::ClimateStatisticsAccumulator uninterrupted_statistics(driver.grid(), 2,
                                                             statistics_start);
  mps::Real uninterrupted_energy = 0.0;
  driver.advance(uninterrupted, parameters.run.end_time_s,
                 [&](const mps::DryHydrostaticState& state,
                     const mps::DryHydrostaticDerived& derived,
                     const mps::DryHydrostaticStepDiagnostics& diagnostics) {
                   uninterrupted_statistics.observe(state, derived);
                   if (state.time_s > statistics_start)
                     uninterrupted_energy +=
                         diagnostics.thermal_energy_contribution_j +
                         diagnostics.rayleigh_drag_energy_contribution_j;
                 });

  auto stopped = driver.initial_state();
  driver.advance(stopped, restart_time);
  auto restarted = mps::unflatten_dry_hydrostatic_state(
      stopped.time_s, stopped.step, mps::flatten_dry_hydrostatic_state(stopped, 2),
      driver.grid().cell_count(), 2);
  mps::ClimateStatisticsAccumulator restarted_statistics(driver.grid(), 2,
                                                         statistics_start);
  mps::Real restarted_energy = 0.0;
  driver.advance(restarted, parameters.run.end_time_s,
                 [&](const mps::DryHydrostaticState& state,
                     const mps::DryHydrostaticDerived& derived,
                     const mps::DryHydrostaticStepDiagnostics& diagnostics) {
                   restarted_statistics.observe(state, derived);
                   if (state.time_s > statistics_start)
                     restarted_energy +=
                         diagnostics.thermal_energy_contribution_j +
                         diagnostics.rayleigh_drag_energy_contribution_j;
                 });

  MPS_CHECK(mps::flatten_dry_hydrostatic_state(uninterrupted, 2) ==
            mps::flatten_dry_hydrostatic_state(restarted, 2));
  MPS_CHECK_EQ(uninterrupted_energy, restarted_energy);
  const auto first = uninterrupted_statistics.rows();
  const auto second = restarted_statistics.rows();
  MPS_CHECK_EQ(first.size(), second.size());
  for (std::size_t index = 0; index < first.size(); ++index) {
    MPS_CHECK_EQ(first[index].mean_temperature_k, second[index].mean_temperature_k);
    MPS_CHECK_EQ(first[index].eddy_kinetic_energy_m2_s2,
                 second[index].eddy_kinetic_energy_m2_s2);
    MPS_CHECK_EQ(first[index].mass_streamfunction_kg_s,
                 second[index].mass_streamfunction_kg_s);
    MPS_CHECK_EQ(first[index].accumulated_time_s, second[index].accumulated_time_s);
  }
}

int main() { return mps::test::run_all(); }
