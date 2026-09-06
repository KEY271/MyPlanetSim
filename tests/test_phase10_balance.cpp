#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <string_view>

#include "myplanetsim/config/experiment_config.hpp"
#include "myplanetsim/dynamics/dry_hydrostatic_driver.hpp"
#include "support/test.hpp"

namespace {

[[nodiscard]] mps::ExperimentConfig semi_implicit_config(const std::string_view name,
                                                         const mps::Real time_step_s,
                                                         const mps::Real duration_s) {
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
      .implicit_weight = 0.5,
      .wave_cfl_threshold = 0.45,
      .maximum_implicit_modes = std::min<mps::Index>(5, config.vertical.levels),
      .nonlinear_iterations = 2,
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

void require_finite_state(const mps::DryHydrostaticState& state) {
  for (const auto value : state.surface_pressure_pa) MPS_CHECK(std::isfinite(value));
  for (const auto value : state.horizontal_momentum_mass_kg_m_s)
    MPS_CHECK(mps::is_finite(value));
  for (const auto value : state.potential_temperature_mass_k_kg_m2)
    MPS_CHECK(std::isfinite(value));
  for (const auto value : state.tracer_mass_kg_m2) MPS_CHECK(std::isfinite(value));
}

}  // namespace

MPS_TEST_CASE("semi-implicit flat and DCMIP terrain rest remain stationary") {
  for (const auto name : {"phase5_isothermal_rest.cfg", "phase6_dcmip_2_0_0.cfg"}) {
    const auto config = semi_implicit_config(name, 1800.0, 3600.0);
    mps::DryHydrostaticDriver driver(config);
    const auto initial = driver.initial_state();
    auto final = initial;
    driver.advance(final, config.run.end_time_s);
    MPS_CHECK(final.surface_pressure_pa == initial.surface_pressure_pa);
    for (std::size_t n = 0; n < final.horizontal_momentum_mass_kg_m_s.size(); ++n)
      MPS_CHECK_EQ(mps::norm(final.horizontal_momentum_mass_kg_m_s[n] -
                             initial.horizontal_momentum_mass_kg_m_s[n]),
                   0.0);
    MPS_CHECK(final.potential_temperature_mass_k_kg_m2 ==
              initial.potential_temperature_mass_k_kg_m2);
    MPS_CHECK(final.tracer_mass_kg_m2 == initial.tracer_mass_kg_m2);
  }
}

MPS_TEST_CASE("semi-implicit mountain wave converges under step refinement") {
  constexpr mps::Real duration_s = 7200.0;
  const auto coarse =
      run(semi_implicit_config("phase6_linear_mountain_wave.cfg", 1800.0, duration_s));
  const auto medium =
      run(semi_implicit_config("phase6_linear_mountain_wave.cfg", 900.0, duration_s));
  const auto fine =
      run(semi_implicit_config("phase6_linear_mountain_wave.cfg", 450.0, duration_s));
  require_finite_state(coarse);
  require_finite_state(medium);
  require_finite_state(fine);
  const auto coarse_difference = normalized_state_difference(coarse, medium);
  const auto fine_difference = normalized_state_difference(medium, fine);
  std::cout << "mountain coarse-medium=" << coarse_difference
            << " medium-fine=" << fine_difference << '\n';
  MPS_CHECK(coarse_difference > 0.0);
  MPS_CHECK(fine_difference < coarse_difference);
}

MPS_TEST_CASE("steady and baroclinic benchmarks stay in the small-step envelope") {
  for (const auto name :
       {"phase5_umjs14_steady.cfg", "phase6_jw06_steady.cfg",
        "phase5_umjs14_baroclinic.cfg", "phase6_jw06_baroclinic.cfg"}) {
    constexpr mps::Real duration_s = 1800.0;
    const auto large = run(semi_implicit_config(name, 1800.0, duration_s));
    const auto refined = run(semi_implicit_config(name, 450.0, duration_s));
    require_finite_state(large);
    require_finite_state(refined);
    const auto difference = normalized_state_difference(large, refined);
    std::cout << name << " large-refined=" << difference << '\n';
    MPS_CHECK(difference < 5.0e-4);
  }
}

int main() { return mps::test::run_all(); }
