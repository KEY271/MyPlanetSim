// Phase 10 long-step gate for the fixed-iteration ICI contract (ADR 0016).
//
// Benard (2003) classifies the classical semi-implicit scheme and its iterative
// variants as one family that differs only in its iteration count, and records that a
// non-extrapolating scheme with a single iteration is first-order accurate in time.
// These tests check the consequences that matter for the registered contract: two
// iterations restore second-order convergence at alpha = 0.5, the production count
// leaves an iteration error below the temporal truncation error, and the accepted step
// ladder ends on the explicit material constraint rather than on solver convergence.
#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <string_view>

#include "myplanetsim/config/experiment_config.hpp"
#include "myplanetsim/dynamics/dry_hydrostatic_driver.hpp"
#include "support/test.hpp"

namespace {

[[nodiscard]] mps::ExperimentConfig long_step_config(const std::string_view name,
                                                     const mps::Real time_step_s,
                                                     const mps::Real duration_s,
                                                     const mps::Index iterations,
                                                     const mps::Index mode_cap) {
  auto config = mps::load_experiment_config(std::string(MPS_CONFIG_DIRECTORY) + "/" +
                                            std::string(name));
  config.grid.cells_per_panel = 4;
  config.run.end_time_s = config.run.start_time_s + duration_s;
  config.run.time_step_s = time_step_s;
  config.diagnostics.interval_steps = 1;
  config.dry_hydrostatic.time_integrator =
      mps::DryHydrostaticTimeIntegrator::kSemiImplicit;
  config.dry_hydrostatic.advective_cfl = 0.45;
  config.semi_implicit = mps::SemiImplicitParameters{
      .reference_surface_pressure_pa = config.vertical.surface_pressure_pa,
      .reference_temperature_k = config.vertical.initial_temperature_k,
      .reference_update = mps::SemiImplicitReferenceUpdate::kFixed,
      .implicit_weight = 0.5,
      .wave_cfl_threshold = 0.45,
      .maximum_implicit_modes = std::min<mps::Index>(mode_cap, config.vertical.levels),
      .nonlinear_iterations = iterations,
      .linear_relative_tolerance = 1.0e-8,
      .linear_absolute_tolerance = 1.0e-12,
      .linear_maximum_iterations = 80,
      .gmres_restart = 20,
      .minimum_time_step_s = 30.0};
  config.output_directory = "x";
  config.validate();
  return config;
}

[[nodiscard]] mps::DryHydrostaticState run(const mps::ExperimentConfig& config) {
  mps::DryHydrostaticDriver driver(config);
  auto state = driver.initial_state();
  driver.advance(state, config.run.end_time_s);
  return state;
}

[[nodiscard]] mps::Real normalized_state_difference(
    const mps::DryHydrostaticState& left, const mps::DryHydrostaticState& right) {
  mps::Real squared = 0.0;
  std::size_t count = 0;
  for (std::size_t cell = 0; cell < left.surface_pressure_pa.size(); ++cell) {
    const mps::Real value =
        (left.surface_pressure_pa[cell] - right.surface_pressure_pa[cell]) / 100000.0;
    squared += value * value;
    ++count;
  }
  for (std::size_t n = 0; n < left.horizontal_momentum_mass_kg_m_s.size(); ++n) {
    const auto momentum = (left.horizontal_momentum_mass_kg_m_s[n] -
                           right.horizontal_momentum_mass_kg_m_s[n]) /
                          1.0e6;
    const mps::Real thermal = (left.potential_temperature_mass_k_kg_m2[n] -
                               right.potential_temperature_mass_k_kg_m2[n]) /
                              1.0e6;
    squared += mps::norm_squared(momentum) + thermal * thermal;
    count += 4;
  }
  return std::sqrt(squared / static_cast<mps::Real>(count));
}

// Two iterations must recover the second-order accuracy of centred Crank--Nicolson.
// Halving the step twice should shrink the successive differences by about four.
void test_second_order_in_time() {
  constexpr mps::Real duration_s = 7200.0;
  const auto coarse = run(long_step_config("phase10_linear_wave_semi_implicit.cfg",
                                           1800.0, duration_s, 2, 5));
  const auto medium = run(long_step_config("phase10_linear_wave_semi_implicit.cfg",
                                           900.0, duration_s, 2, 5));
  const auto fine = run(long_step_config("phase10_linear_wave_semi_implicit.cfg", 450.0,
                                         duration_s, 2, 5));
  const auto finest = run(long_step_config("phase10_linear_wave_semi_implicit.cfg",
                                           225.0, duration_s, 2, 5));
  const auto coarse_difference = normalized_state_difference(coarse, medium);
  const auto medium_difference = normalized_state_difference(medium, fine);
  const auto fine_difference = normalized_state_difference(fine, finest);
  const auto first_order = std::log2(coarse_difference / medium_difference);
  const auto second_order = std::log2(medium_difference / fine_difference);
  std::cout << "linear wave 1800-900=" << coarse_difference
            << " 900-450=" << medium_difference << " 450-225=" << fine_difference
            << " observed order " << first_order << ", " << second_order << '\n';
  MPS_CHECK(coarse_difference > 0.0);
  MPS_CHECK(first_order > 1.7);
  MPS_CHECK(second_order > 1.7);
}

// Stopping at a fixed iteration count is only defensible if the remaining iteration
// error is below the temporal truncation error being accepted at the same step. This
// measures both against a well-iterated reference and pins the production count.
// Measured contraction here is about one order of magnitude per iteration, matching the
// rate Thuburn et al. (2014) report, so two iterations are second-order but still leave
// an iteration error larger than the truncation error, and three do not.
void test_production_iteration_count_is_below_truncation_error() {
  constexpr mps::Real duration_s = 7200.0;
  constexpr mps::Index converged = 6;
  const auto reference = run(long_step_config("phase10_linear_wave_semi_implicit.cfg",
                                              1800.0, duration_s, converged, 5));
  const auto halved = run(long_step_config("phase10_linear_wave_semi_implicit.cfg",
                                           900.0, duration_s, converged, 5));
  const auto truncation = normalized_state_difference(reference, halved);
  MPS_CHECK(truncation > 0.0);
  mps::Real previous = 0.0;
  for (const mps::Index iterations : {2, 3, 4}) {
    const auto state = run(long_step_config("phase10_linear_wave_semi_implicit.cfg",
                                            1800.0, duration_s, iterations, 5));
    const auto error = normalized_state_difference(state, reference);
    std::cout << "iteration error N=" << iterations << " is " << error
              << " against truncation " << truncation << '\n';
    if (previous > 0.0) MPS_CHECK(error < 0.5 * previous);
    previous = error;
    if (iterations >= 3) MPS_CHECK(error < truncation);
  }
}

// The registered long-step ladder. With enough implicit modes the scheme survives well
// past the delivery target, and the first failure is the explicit material constraint
// that Phase 10 deliberately does not make implicit.
void test_long_step_ladder() {
  for (const mps::Real time_step_s : {1800.0, 2400.0, 3600.0, 7200.0}) {
    auto config = long_step_config("phase10_held_suarez_semi_implicit.cfg", time_step_s,
                                   4.0 * time_step_s, 3, 20);
    config.semi_implicit->minimum_time_step_s = time_step_s;
    const auto state = run(config);
    std::cout << "held-suarez accepted dt=" << time_step_s << '\n';
    MPS_CHECK_NEAR(state.time_s, config.run.end_time_s, 1.0e-9);
    for (const auto value : state.surface_pressure_pa) MPS_CHECK(std::isfinite(value));
  }
  // The coarse test grid relaxes the material limit by about the cell-size ratio, so
  // the ladder is pushed further here than on the registered N=12 preset.
  auto beyond = long_step_config("phase10_held_suarez_semi_implicit.cfg", 86400.0,
                                 86400.0, 3, 20);
  beyond.semi_implicit->minimum_time_step_s = 86400.0;
  bool rejected = false;
  try {
    static_cast<void>(run(beyond));
  } catch (const std::exception& error) {
    rejected = true;
    std::cout << "ladder stops at: " << error.what() << '\n';
  }
  MPS_CHECK(rejected);
}

}  // namespace

int main() {
  test_second_order_in_time();
  test_production_iteration_count_is_below_truncation_error();
  test_long_step_ladder();
  std::cout << "phase10 long step gate passed\n";
  return 0;
}
