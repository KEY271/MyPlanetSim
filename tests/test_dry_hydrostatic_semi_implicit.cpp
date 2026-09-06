#include <algorithm>
#include <cmath>
#include <vector>

#include "myplanetsim/dynamics/dry_hydrostatic_driver.hpp"
#include "myplanetsim/dynamics/dry_hydrostatic_semi_implicit.hpp"
#include "support/test.hpp"

namespace {

constexpr mps::PlanetParameters planet{6371220.0, 7.29212e-5, 9.80616,
                                       287.0,     1004.0,     100000.0};

[[nodiscard]] mps::AtmosphericHybridCoordinate coordinate() {
  return mps::AtmosphericHybridCoordinate(
      {.a_half_pa = {1000.0, 750.0, 500.0, 250.0, 0.0},
       .b_half = {0.0, 0.25, 0.5, 0.75, 1.0}},
      80000.0, 120000.0, 100.0);
}

constexpr mps::SemiImplicitParameters parameters{
    .reference_surface_pressure_pa = 100000.0,
    .reference_temperature_k = 280.0,
    .implicit_weight = 0.5,
    .wave_cfl_threshold = 0.45,
    .maximum_implicit_modes = 4,
    .nonlinear_relative_tolerance = 1.0e-8,
    .nonlinear_maximum_iterations = 4,
    .linear_relative_tolerance = 1.0e-10,
    .linear_absolute_tolerance = 1.0e-12,
    .linear_maximum_iterations = 80,
    .gmres_restart = 20,
    .minimum_time_step_s = 60.0};

[[nodiscard]] mps::ExperimentConfig long_step_config() {
  const auto coefficients = mps::uniform_sigma_coefficients(1000.0, 20);
  mps::ExperimentConfig result{
      .kind = mps::ExperimentKind::kDryHydrostatic,
      .planet = planet,
      .run = {0.0, 1800.0, 1800.0, 0},
      .grid = {12},
      .vertical = {.levels = 20,
                   .a_half_pa = coefficients.a_half_pa,
                   .b_half = coefficients.b_half,
                   .surface_pressure_pa = 100000.0,
                   .minimum_surface_pressure_pa = 80000.0,
                   .maximum_surface_pressure_pa = 120000.0,
                   .minimum_pressure_thickness_pa = 100.0,
                   .initial_temperature_k = 288.0,
                   .initial_potential_temperature_k = 300.0,
                   .temperature_floor_k = 100.0,
                   .transport_scheme = mps::VerticalTransportScheme::kLinear,
                   .limiter = mps::VerticalLimiterKind::kMinmod,
                   .cfl = 0.5},
      .dry_hydrostatic = {.test_case = mps::DryHydrostaticTestCase::kLinearWave,
                          .reconstruction = mps::ReconstructionKind::kLinear,
                          .limiter = mps::LimiterKind::kBarthJespersen,
                          .cfl = 0.45,
                          .diffusion_kind = mps::DiffusionKind::kNone,
                          .diffusion_coefficient = 0.0,
                          .time_integrator =
                              mps::DryHydrostaticTimeIntegrator::kSemiImplicit,
                          .advective_cfl = 0.45},
      .semi_implicit = parameters,
      .diagnostics = {1},
      .output_directory = "x"};
  result.semi_implicit->reference_temperature_k = 288.0;
  result.validate();
  return result;
}

}  // namespace

MPS_TEST_CASE("external-mode Helmholtz solve inverts the coupled correction") {
  const auto c = coordinate();
  const auto reference =
      mps::make_dry_hydrostatic_reference_column(c, planet, parameters);
  const auto modes = mps::make_dry_hydrostatic_external_mode(reference, planet);
  const auto fast = mps::make_dry_hydrostatic_fast_operator(c, planet, reference);
  const auto external =
      mps::make_dry_hydrostatic_external_mode_operator(reference, modes);
  const mps::CubedSphereGrid grid(4, planet.radius_m);
  const auto cells = grid.cell_count();
  const auto volume = cells * fast.levels;
  mps::DryHydrostaticFastPerturbation right_hand_side{
      .surface_pressure_pa = std::vector<mps::Real>(cells),
      .horizontal_momentum_mass_kg_m_s = std::vector<mps::Vec3>(volume),
      .potential_temperature_mass_k_kg_m2 = std::vector<mps::Real>(volume)};
  for (std::size_t cell = 0; cell < cells; ++cell) {
    const auto center = grid.cells()[cell].center;
    right_hand_side.surface_pressure_pa[cell] = 10.0 * (center.x - 0.2 * center.y);
    const auto tangent = mps::project_tangent({-center.y, center.x, 0.1}, center);
    for (std::size_t level = 0; level < fast.levels; ++level) {
      const auto n = mps::dry_hydrostatic_offset(cell, level, fast.levels);
      right_hand_side.horizontal_momentum_mass_kg_m_s[n] =
          (1.0 + static_cast<double>(level)) * tangent;
      right_hand_side.potential_temperature_mass_k_kg_m2[n] =
          0.5 * (1.0 + static_cast<double>(level)) * center.z;
    }
  }
  mps::DryHydrostaticSemiImplicitWorkspace workspace;
  mps::DryHydrostaticFastPerturbation correction;
  const auto solve = mps::solve_dry_hydrostatic_external_mode_correction(
      grid, planet, fast, external, 900.0, right_hand_side,
      {.restart = 20,
       .maximum_iterations = 80,
       .relative_tolerance = 1.0e-10,
       .absolute_tolerance = 1.0e-12},
      correction, workspace);
  MPS_CHECK(solve.linear.converged());
  MPS_CHECK(solve.equation_residual_norm < 2.0e-9);
  MPS_CHECK(solve.linear.iterations < 80);
}

MPS_TEST_CASE("external-mode operator preserves fast scalar integrals") {
  const auto c = coordinate();
  const auto reference =
      mps::make_dry_hydrostatic_reference_column(c, planet, parameters);
  const auto modes = mps::make_dry_hydrostatic_external_mode(reference, planet);
  const auto fast = mps::make_dry_hydrostatic_fast_operator(c, planet, reference);
  const auto external =
      mps::make_dry_hydrostatic_external_mode_operator(reference, modes);
  const mps::CubedSphereGrid grid(3, planet.radius_m);
  mps::DryHydrostaticFastPerturbation delta{
      .surface_pressure_pa = std::vector<mps::Real>(grid.cell_count()),
      .horizontal_momentum_mass_kg_m_s =
          std::vector<mps::Vec3>(grid.cell_count() * fast.levels),
      .potential_temperature_mass_k_kg_m2 =
          std::vector<mps::Real>(grid.cell_count() * fast.levels)};
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
    const auto center = grid.cells()[cell].center;
    delta.surface_pressure_pa[cell] = 20.0 * center.x;
    for (std::size_t level = 0; level < fast.levels; ++level) {
      const auto n = mps::dry_hydrostatic_offset(cell, level, fast.levels);
      delta.horizontal_momentum_mass_kg_m_s[n] =
          fast.reference_air_mass_kg_m2[level] *
          mps::project_tangent({-center.y, center.x, 0.3}, center);
    }
  }
  mps::DryHydrostaticFastOperatorWorkspace workspace;
  mps::DryHydrostaticFastTendency tendency;
  mps::apply_dry_hydrostatic_external_mode_operator(grid, planet, fast, external, delta,
                                                    tendency, workspace);
  mps::Real pressure = 0.0;
  mps::Real pressure_scale = 0.0;
  mps::Real theta = 0.0;
  mps::Real theta_scale = 0.0;
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
    const auto area = grid.cells()[cell].area_m2;
    pressure += area * tendency.surface_pressure_pa_s[cell];
    pressure_scale += area * std::abs(tendency.surface_pressure_pa_s[cell]);
    for (std::size_t level = 0; level < fast.levels; ++level) {
      const auto n = mps::dry_hydrostatic_offset(cell, level, fast.levels);
      theta += area * tendency.tendency.potential_temperature_mass[n];
      theta_scale += area * std::abs(tendency.tendency.potential_temperature_mass[n]);
    }
  }
  MPS_CHECK_NEAR(pressure, 0.0, 2.0e-13 * std::max(1.0, pressure_scale));
  MPS_CHECK_NEAR(theta, 0.0, 2.0e-13 * std::max(1.0, theta_scale));
}

MPS_TEST_CASE("dry linear wave accepts one 1800-second Crank-Nicolson step") {
  const auto config = long_step_config();
  mps::DryHydrostaticDriver driver(config);
  auto state = driver.initial_state();
  const auto initial_rhs = driver.rhs(state);
  MPS_CHECK(initial_rhs.horizontal_fast_wave_stable_time_step_s <
            config.run.time_step_s);
  const auto initial_pressure = state.surface_pressure_pa;
  const auto scalar_integrals = [&](const mps::DryHydrostaticState& value) {
    mps::Real pressure = 0.0;
    mps::Real potential_temperature_mass = 0.0;
    for (std::size_t cell = 0; cell < driver.grid().cell_count(); ++cell) {
      const auto area = driver.grid().cells()[cell].area_m2;
      pressure += area * value.surface_pressure_pa[cell];
      for (std::size_t level = 0; level < 20; ++level) {
        potential_temperature_mass +=
            area * value.potential_temperature_mass_k_kg_m2[mps::dry_hydrostatic_offset(
                       cell, level, 20)];
      }
    }
    return std::pair{pressure, potential_temperature_mass};
  };
  const auto initial_integrals = scalar_integrals(state);
  driver.advance(state, config.run.end_time_s);
  MPS_CHECK_EQ(state.step, 1U);
  MPS_CHECK_NEAR(state.time_s, 1800.0, 0.0);
  mps::Real maximum_pressure_change = 0.0;
  for (std::size_t cell = 0; cell < state.surface_pressure_pa.size(); ++cell) {
    MPS_CHECK(std::isfinite(state.surface_pressure_pa[cell]));
    maximum_pressure_change =
        std::max(maximum_pressure_change,
                 std::abs(state.surface_pressure_pa[cell] - initial_pressure[cell]));
  }
  MPS_CHECK(maximum_pressure_change > 0.0);
  const auto final_integrals = scalar_integrals(state);
  MPS_CHECK_NEAR(final_integrals.first, initial_integrals.first,
                 2.0e-12 * std::abs(initial_integrals.first));
  MPS_CHECK_NEAR(final_integrals.second, initial_integrals.second,
                 2.0e-12 * std::abs(initial_integrals.second));
}

int main() { return mps::test::run_all(); }
