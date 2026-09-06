#include <algorithm>
#include <cmath>
#include <vector>

#include "myplanetsim/dynamics/dry_hydrostatic_fast_operator.hpp"
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
    .nonlinear_iterations = 2,
    .linear_relative_tolerance = 1.0e-8,
    .linear_absolute_tolerance = 1.0e-12,
    .linear_maximum_iterations = 40,
    .gmres_restart = 20,
    .minimum_time_step_s = 60.0};

[[nodiscard]] mps::DryHydrostaticFastPerturbation perturbation(
    const mps::CubedSphereGrid& grid,
    const mps::DryHydrostaticReferenceColumn& reference, const double shift) {
  const auto cells = grid.cell_count();
  const auto levels = reference.geometry.air_mass_kg_m2.size();
  mps::DryHydrostaticFastPerturbation result{
      .surface_pressure_pa = std::vector<mps::Real>(cells),
      .horizontal_momentum_mass_kg_m_s = std::vector<mps::Vec3>(cells * levels),
      .potential_temperature_mass_k_kg_m2 = std::vector<mps::Real>(cells * levels)};
  for (std::size_t cell = 0; cell < cells; ++cell) {
    const auto center = grid.cells()[cell].center;
    result.surface_pressure_pa[cell] =
        20.0 * (center.x + 0.3 * center.y + shift * center.z);
    const auto direction =
        mps::project_tangent({-center.y + shift, center.x, 0.2}, center);
    for (std::size_t level = 0; level < levels; ++level) {
      const auto n = mps::dry_hydrostatic_offset(cell, level, levels);
      result.horizontal_momentum_mass_kg_m_s[n] =
          (3.0 + static_cast<double>(level)) *
          reference.geometry.air_mass_kg_m2[level] * direction;
      result.potential_temperature_mass_k_kg_m2[n] =
          1.0e-4 * reference.potential_temperature_mass_k_kg_m2[level] *
          (center.z + shift * center.x);
    }
  }
  return result;
}

[[nodiscard]] mps::DryHydrostaticFastPerturbation linear_combination(
    const mps::DryHydrostaticFastPerturbation& a, const double a_scale,
    const mps::DryHydrostaticFastPerturbation& b, const double b_scale) {
  auto result = a;
  for (std::size_t cell = 0; cell < result.surface_pressure_pa.size(); ++cell)
    result.surface_pressure_pa[cell] =
        a_scale * a.surface_pressure_pa[cell] + b_scale * b.surface_pressure_pa[cell];
  for (std::size_t n = 0; n < result.horizontal_momentum_mass_kg_m_s.size(); ++n) {
    result.horizontal_momentum_mass_kg_m_s[n] =
        a_scale * a.horizontal_momentum_mass_kg_m_s[n] +
        b_scale * b.horizontal_momentum_mass_kg_m_s[n];
    result.potential_temperature_mass_k_kg_m2[n] =
        a_scale * a.potential_temperature_mass_k_kg_m2[n] +
        b_scale * b.potential_temperature_mass_k_kg_m2[n];
  }
  return result;
}

void check_linear_combination(const mps::DryHydrostaticFastTendency& actual,
                              const mps::DryHydrostaticFastTendency& a,
                              const double a_scale,
                              const mps::DryHydrostaticFastTendency& b,
                              const double b_scale) {
  for (std::size_t cell = 0; cell < actual.surface_pressure_pa_s.size(); ++cell) {
    const double expected = a_scale * a.surface_pressure_pa_s[cell] +
                            b_scale * b.surface_pressure_pa_s[cell];
    MPS_CHECK_NEAR(actual.surface_pressure_pa_s[cell], expected,
                   2.0e-12 * std::max(1.0, std::abs(expected)));
  }
  for (std::size_t n = 0; n < actual.tendency.air_mass.size(); ++n) {
    const double expected_mass =
        a_scale * a.tendency.air_mass[n] + b_scale * b.tendency.air_mass[n];
    const double expected_theta = a_scale * a.tendency.potential_temperature_mass[n] +
                                  b_scale * b.tendency.potential_temperature_mass[n];
    const auto expected_momentum =
        a_scale * a.tendency.momentum[n] + b_scale * b.tendency.momentum[n];
    MPS_CHECK_NEAR(actual.tendency.air_mass[n], expected_mass,
                   2.0e-12 * std::max(1.0, std::abs(expected_mass)));
    MPS_CHECK_NEAR(actual.tendency.potential_temperature_mass[n], expected_theta,
                   2.0e-12 * std::max(1.0, std::abs(expected_theta)));
    MPS_CHECK_NEAR(mps::norm(actual.tendency.momentum[n] - expected_momentum), 0.0,
                   2.0e-12 * std::max(1.0, mps::norm(expected_momentum)));
    MPS_CHECK_EQ(actual.tendency.tracer_mass[n], 0.0);
  }
}

}  // namespace

MPS_TEST_CASE("reference-linear fast operator is linear and its rest state is zero") {
  const auto c = coordinate();
  const auto reference =
      mps::make_dry_hydrostatic_reference_column(c, planet, parameters);
  const auto op = mps::make_dry_hydrostatic_fast_operator(c, planet, reference);
  const mps::CubedSphereGrid grid(2, planet.radius_m);
  const auto a = perturbation(grid, reference, 0.2);
  const auto b = perturbation(grid, reference, -0.7);
  const auto combined = linear_combination(a, 1.7, b, -0.4);
  mps::DryHydrostaticFastOperatorWorkspace workspace;
  mps::DryHydrostaticFastTendency la;
  mps::DryHydrostaticFastTendency lb;
  mps::DryHydrostaticFastTendency actual;
  mps::DryHydrostaticFastTendency zero;
  mps::DryHydrostaticFastPerturbation zero_perturbation{
      .surface_pressure_pa = std::vector<mps::Real>(grid.cell_count()),
      .horizontal_momentum_mass_kg_m_s =
          std::vector<mps::Vec3>(grid.cell_count() * op.levels),
      .potential_temperature_mass_k_kg_m2 =
          std::vector<mps::Real>(grid.cell_count() * op.levels)};
  mps::apply_dry_hydrostatic_fast_operator(grid, planet, op, a, la, workspace);
  mps::apply_dry_hydrostatic_fast_operator(grid, planet, op, b, lb, workspace);
  mps::apply_dry_hydrostatic_fast_operator(grid, planet, op, combined, actual,
                                           workspace);
  mps::apply_dry_hydrostatic_fast_operator(grid, planet, op, zero_perturbation, zero,
                                           workspace);
  check_linear_combination(actual, la, 1.7, lb, -0.4);
  for (const auto value : zero.surface_pressure_pa_s) MPS_CHECK_EQ(value, 0.0);
  for (std::size_t n = 0; n < zero.tendency.air_mass.size(); ++n) {
    MPS_CHECK_EQ(zero.tendency.air_mass[n], 0.0);
    MPS_CHECK_EQ(mps::norm(zero.tendency.momentum[n]), 0.0);
    MPS_CHECK_EQ(zero.tendency.potential_temperature_mass[n], 0.0);
  }
}

MPS_TEST_CASE("fast mass and potential-temperature fluxes conserve global integrals") {
  const auto c = coordinate();
  const auto reference =
      mps::make_dry_hydrostatic_reference_column(c, planet, parameters);
  const auto op = mps::make_dry_hydrostatic_fast_operator(c, planet, reference);
  const mps::CubedSphereGrid grid(3, planet.radius_m);
  const auto delta = perturbation(grid, reference, 0.35);
  mps::DryHydrostaticFastOperatorWorkspace workspace;
  mps::DryHydrostaticFastTendency tendency;
  mps::apply_dry_hydrostatic_fast_operator(grid, planet, op, delta, tendency,
                                           workspace);

  mps::Real pressure_integral = 0.0;
  mps::Real pressure_absolute = 0.0;
  mps::Real theta_integral = 0.0;
  mps::Real theta_absolute = 0.0;
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
    const auto area = grid.cells()[cell].area_m2;
    pressure_integral += area * tendency.surface_pressure_pa_s[cell];
    pressure_absolute += area * std::abs(tendency.surface_pressure_pa_s[cell]);
    for (std::size_t level = 0; level < op.levels; ++level) {
      const auto n = mps::dry_hydrostatic_offset(cell, level, op.levels);
      theta_integral += area * tendency.tendency.potential_temperature_mass[n];
      theta_absolute +=
          area * std::abs(tendency.tendency.potential_temperature_mass[n]);
    }
  }
  MPS_CHECK_NEAR(pressure_integral, 0.0, 2.0e-13 * std::max(1.0, pressure_absolute));
  MPS_CHECK_NEAR(theta_integral, 0.0, 2.0e-13 * std::max(1.0, theta_absolute));
}

MPS_TEST_CASE("fast add-subtract split reconstructs the original tendency") {
  const auto c = coordinate();
  const auto reference =
      mps::make_dry_hydrostatic_reference_column(c, planet, parameters);
  const auto op = mps::make_dry_hydrostatic_fast_operator(c, planet, reference);
  const mps::CubedSphereGrid grid(2, planet.radius_m);
  mps::DryHydrostaticFastOperatorWorkspace workspace;
  mps::DryHydrostaticFastTendency full;
  mps::DryHydrostaticFastTendency linear;
  mps::apply_dry_hydrostatic_fast_operator(
      grid, planet, op, perturbation(grid, reference, 0.6), full, workspace);
  mps::apply_dry_hydrostatic_fast_operator(
      grid, planet, op, perturbation(grid, reference, -0.25), linear, workspace);
  const auto remainder = mps::subtract_dry_hydrostatic_fast_tendency(full, linear);
  const auto reconstructed = mps::add_dry_hydrostatic_fast_tendency(linear, remainder);
  check_linear_combination(reconstructed, full, 1.0, full, 0.0);
}

MPS_TEST_CASE("state subtraction maps the horizontal reference state to zero") {
  const auto c = coordinate();
  const auto reference =
      mps::make_dry_hydrostatic_reference_column(c, planet, parameters);
  constexpr std::size_t cells = 3;
  mps::DryHydrostaticState state{
      .surface_pressure_pa =
          std::vector<mps::Real>(cells, reference.surface_pressure_pa),
      .horizontal_momentum_mass_kg_m_s = std::vector<mps::Vec3>(cells * c.levels()),
      .potential_temperature_mass_k_kg_m2 = std::vector<mps::Real>(cells * c.levels()),
      .tracer_mass_kg_m2 = std::vector<mps::Real>(cells * c.levels()),
      .surface_temperature_k = {}};
  for (std::size_t cell = 0; cell < cells; ++cell)
    std::copy(reference.potential_temperature_mass_k_kg_m2.begin(),
              reference.potential_temperature_mass_k_kg_m2.end(),
              state.potential_temperature_mass_k_kg_m2.begin() +
                  static_cast<std::ptrdiff_t>(cell * c.levels()));
  const auto delta = mps::make_dry_hydrostatic_fast_perturbation(state, reference);
  MPS_CHECK(std::ranges::all_of(delta.surface_pressure_pa,
                                [](const auto value) { return value == 0.0; }));
  MPS_CHECK(std::ranges::all_of(delta.potential_temperature_mass_k_kg_m2,
                                [](const auto value) { return value == 0.0; }));
}

int main() { return mps::test::run_all(); }
