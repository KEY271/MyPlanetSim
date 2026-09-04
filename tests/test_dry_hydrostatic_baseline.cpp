#include <cmath>
#include <sstream>

#include "myplanetsim/diagnostics/dry_hydrostatic_diagnostics.hpp"
#include "myplanetsim/dynamics/dry_hydrostatic_driver.hpp"
#include "myplanetsim/io/checkpoint.hpp"
#include "myplanetsim/io/run_metadata.hpp"
#include "support/test.hpp"

namespace {

// The registered unforced, flat-surface baseline is the Phase 5 UMJS14 steady preset
// at CI scale. Phase 7 compares its forced runs against these values, so the
// configuration, window, and tolerances below are fixed before any physics is added.
constexpr double kEndTimeS = 600.0;
constexpr double kRestartTimeS = 300.0;
// Observed relative drift over the window is at rounding for the masses, 5.4e-7 for
// total energy, and 1.3e-5 for absolute axial angular momentum. The
// Rusanov/least-squares reference scheme does not conserve energy or angular momentum
// exactly, so these bound the drift instead of asserting conservation.
constexpr double kMassDriftTolerance = 1.0e-13;
constexpr double kEnergyDriftTolerance = 1.0e-6;
constexpr double kAngularMomentumDriftTolerance = 2.0e-5;

[[nodiscard]] mps::ExperimentConfig baseline_config() {
  mps::ExperimentConfig config;
  config.kind = mps::ExperimentKind::kDryHydrostatic;
  config.planet = {6371220.0, 7.29212e-5, 9.80616, 287.0, 1004.0, 100000.0};
  config.run = {0.0, kEndTimeS, 10.0, 0};
  config.grid = {4};
  config.vertical.levels = 8;
  config.vertical.a_half_pa = {1000, 875, 750, 625, 500, 375, 250, 125, 0};
  config.vertical.b_half = {0, .125, .25, .375, .5, .625, .75, .875, 1};
  config.vertical.surface_pressure_pa = 100000.0;
  config.vertical.minimum_surface_pressure_pa = 80000.0;
  config.vertical.maximum_surface_pressure_pa = 120000.0;
  config.vertical.minimum_pressure_thickness_pa = 100.0;
  config.vertical.initial_temperature_k = 288.0;
  config.vertical.initial_potential_temperature_k = 300.0;
  config.vertical.temperature_floor_k = 100.0;
  config.vertical.transport_scheme = mps::VerticalTransportScheme::kLinear;
  config.vertical.limiter = mps::VerticalLimiterKind::kMinmod;
  config.vertical.cfl = 0.5;
  config.dry_hydrostatic.test_case = mps::DryHydrostaticTestCase::kUmjs14Steady;
  config.output_directory = "output";
  config.validate();
  return config;
}

[[nodiscard]] mps::DryHydrostaticDiagnostics budgets(
    const mps::DryHydrostaticDriver& driver, const mps::ExperimentConfig& config,
    const mps::DryHydrostaticState& state) {
  return diagnose_dry_hydrostatic_budgets(driver.grid(), state, driver.diagnose(state),
                                          config.planet);
}

[[nodiscard]] double relative_drift(const double initial, const double final_value) {
  return std::abs(final_value - initial) / std::max(std::abs(initial), 1.0);
}

// Round trips through the shipped checkpoint payload so the restart baseline covers the
// serialized state rather than an in-memory copy.
[[nodiscard]] mps::DryHydrostaticState round_trip(
    const mps::ExperimentConfig& config, const mps::DryHydrostaticState& state) {
  const auto fingerprint = mps::config_fingerprint(config);
  const auto levels = static_cast<std::size_t>(config.vertical.levels);
  const auto flat = mps::flatten_dry_hydrostatic_state(state, levels);
  std::stringstream stream;
  mps::write_checkpoint(
      stream, {.time_s = state.time_s,
               .step = state.step,
               .state = flat,
               .config_fingerprint = fingerprint,
               .layout_id = std::string(mps::kDryHydrostaticCheckpointLayout)});
  const auto checkpoint = mps::read_checkpoint(
      stream, fingerprint, mps::kDryHydrostaticCheckpointLayout, flat.size());
  return mps::unflatten_dry_hydrostatic_state(checkpoint.time_s, checkpoint.step,
                                              checkpoint.state,
                                              state.surface_pressure_pa.size(), levels);
}

}  // namespace

MPS_TEST_CASE("the unforced flat baseline keeps mass, energy, and angular momentum") {
  const auto config = baseline_config();
  const mps::DryHydrostaticDriver driver(config);
  auto state = driver.initial_state();
  const auto initial = budgets(driver, config, state);
  driver.advance(state, kEndTimeS);
  const auto final_budgets = budgets(driver, config, state);

  MPS_CHECK_EQ(state.time_s, kEndTimeS);
  MPS_CHECK(state.step > 0);
  MPS_CHECK(relative_drift(initial.dry_mass_kg, final_budgets.dry_mass_kg) <
            kMassDriftTolerance);
  MPS_CHECK(relative_drift(initial.potential_temperature_mass_k_kg,
                           final_budgets.potential_temperature_mass_k_kg) <
            kMassDriftTolerance);
  MPS_CHECK(relative_drift(initial.tracer_mass_kg, final_budgets.tracer_mass_kg) <
            kMassDriftTolerance);
  MPS_CHECK(relative_drift(initial.total_energy_j, final_budgets.total_energy_j) <
            kEnergyDriftTolerance);
  MPS_CHECK(relative_drift(initial.axial_angular_momentum_kg_m2_s,
                           final_budgets.axial_angular_momentum_kg_m2_s) <
            kAngularMomentumDriftTolerance);
}

MPS_TEST_CASE("a mid-window checkpoint restart reproduces the baseline exactly") {
  const auto config = baseline_config();
  const mps::DryHydrostaticDriver driver(config);
  auto uninterrupted = driver.initial_state();
  driver.advance(uninterrupted, kEndTimeS);

  auto restarted = driver.initial_state();
  driver.advance(restarted, kRestartTimeS);
  MPS_CHECK_EQ(restarted.time_s, kRestartTimeS);
  restarted = round_trip(config, restarted);
  driver.advance(restarted, kEndTimeS);

  const auto levels = static_cast<std::size_t>(config.vertical.levels);
  const auto expected = mps::flatten_dry_hydrostatic_state(uninterrupted, levels);
  const auto actual = mps::flatten_dry_hydrostatic_state(restarted, levels);
  MPS_CHECK_EQ(restarted.step, uninterrupted.step);
  MPS_CHECK_EQ(restarted.time_s, uninterrupted.time_s);
  MPS_CHECK_EQ(actual.size(), expected.size());
  for (std::size_t n = 0; n < actual.size(); ++n) {
    MPS_CHECK_EQ(actual[n], expected[n]);
  }
}

int main() { return mps::test::run_all(); }
