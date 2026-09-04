#include <cmath>
#include <numbers>

#include "myplanetsim/diagnostics/dry_hydrostatic_diagnostics.hpp"
#include "myplanetsim/dynamics/dry_hydrostatic_driver.hpp"
#include "myplanetsim/physics/held_suarez.hpp"
#include "support/test.hpp"

namespace {

[[nodiscard]] mps::ExperimentConfig config() {
  return {.kind = mps::ExperimentKind::kDryHydrostatic,
          .planet = {6371220, 7.29212e-5, 9.80616, 287, 1004, 100000},
          .run = {0, 600, 10, 0},
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

[[nodiscard]] mps::AtmosphericHybridCoordinate coordinate(
    const mps::ExperimentConfig& parameters) {
  return {{parameters.vertical.a_half_pa, parameters.vertical.b_half},
          parameters.vertical.minimum_surface_pressure_pa,
          parameters.vertical.maximum_surface_pressure_pa,
          parameters.vertical.minimum_pressure_thickness_pa};
}

}  // namespace

MPS_TEST_CASE("Held-Suarez rates follow the fixed analytic profile") {
  const auto planet = config().planet;
  const auto equatorial_surface = mps::held_suarez_rates(0.0, 100000, 100000, planet);
  MPS_CHECK_NEAR(equatorial_surface.equilibrium_temperature_k, 315.0, 1e-13);
  MPS_CHECK_NEAR(equatorial_surface.temperature_relaxation_rate_s_1,
                 1.0 / (4.0 * 86400.0), 1e-18);
  MPS_CHECK_NEAR(equatorial_surface.rayleigh_drag_rate_s_1, 1.0 / 86400.0, 1e-18);

  const auto polar_surface = mps::held_suarez_rates(std::numbers::pi_v<mps::Real> / 2.0,
                                                    100000, 100000, planet);
  MPS_CHECK_NEAR(polar_surface.equilibrium_temperature_k, 255.0, 1e-13);
  MPS_CHECK_NEAR(polar_surface.temperature_relaxation_rate_s_1, 1.0 / (40.0 * 86400.0),
                 1e-18);

  const auto boundary = mps::held_suarez_rates(0.0, 70000, 100000, planet);
  MPS_CHECK_NEAR(boundary.rayleigh_drag_rate_s_1, 0.0, 0.0);
  MPS_CHECK_NEAR(boundary.temperature_relaxation_rate_s_1, 1.0 / (40.0 * 86400.0),
                 1e-18);
  const auto aloft = mps::held_suarez_rates(0.0, 1000, 100000, planet);
  MPS_CHECK_NEAR(aloft.equilibrium_temperature_k, 200.0, 0.0);
  MPS_CHECK_NEAR(aloft.rayleigh_drag_rate_s_1, 0.0, 0.0);
}

MPS_TEST_CASE("Held-Suarez tendency preserves mass shapes and tangent momentum") {
  const auto parameters = config();
  const mps::DryHydrostaticDriver driver(parameters);
  const auto state = driver.initial_state();
  auto derived = driver.diagnose(state);
  for (std::size_t cell = 0; cell < derived.cells; ++cell) {
    for (std::size_t level = 0; level < derived.levels; ++level) {
      const auto offset = mps::dry_hydrostatic_offset(cell, level, derived.levels);
      derived.velocity_m_s[offset] =
          20.0 * mps::cross(mps::Vec3{0, 0, 1}, driver.grid().cells()[cell].center);
    }
  }
  const auto forcing =
      mps::held_suarez_tendency(driver.grid(), coordinate(parameters), derived,
                                state.surface_pressure_pa, parameters.planet);
  MPS_CHECK_EQ(forcing.potential_temperature_mass_k_kg_m2_s.size(),
               derived.cells * derived.levels);
  MPS_CHECK(forcing.diagnostics.rayleigh_drag_work_w <= 0.0);
  for (std::size_t cell = 0; cell < derived.cells; ++cell) {
    for (std::size_t level = 0; level < derived.levels; ++level) {
      const auto offset = mps::dry_hydrostatic_offset(cell, level, derived.levels);
      MPS_CHECK_NEAR(mps::dot(forcing.horizontal_momentum_mass_kg_m_s2[offset],
                              driver.grid().cells()[cell].center),
                     0.0, 1e-12);
      if (derived.pressure_pa[offset] / state.surface_pressure_pa[cell] <= 0.7)
        MPS_CHECK_EQ(mps::norm(forcing.horizontal_momentum_mass_kg_m_s2[offset]), 0.0);
    }
  }
}

MPS_TEST_CASE("Held-Suarez equilibrium at rest has exact zero sources") {
  const auto parameters = config();
  const mps::DryHydrostaticDriver driver(parameters);
  const auto state = driver.initial_state();
  auto derived = driver.diagnose(state);
  for (std::size_t cell = 0; cell < derived.cells; ++cell) {
    const auto latitude = std::asin(driver.grid().cells()[cell].center.z);
    for (std::size_t level = 0; level < derived.levels; ++level) {
      const auto offset = mps::dry_hydrostatic_offset(cell, level, derived.levels);
      derived.temperature_k[offset] =
          mps::held_suarez_rates(latitude, derived.pressure_pa[offset], 100000,
                                 parameters.planet)
              .equilibrium_temperature_k;
      derived.velocity_m_s[offset] = {};
    }
  }
  const auto forcing =
      mps::held_suarez_tendency(driver.grid(), coordinate(parameters), derived,
                                state.surface_pressure_pa, parameters.planet);
  for (const auto value : forcing.potential_temperature_mass_k_kg_m2_s)
    MPS_CHECK_EQ(value, 0.0);
  for (const auto value : forcing.horizontal_momentum_mass_kg_m_s2)
    MPS_CHECK_EQ(mps::norm(value), 0.0);
  MPS_CHECK_EQ(forcing.diagnostics.thermal_energy_rate_w, 0.0);
  MPS_CHECK_EQ(forcing.diagnostics.rayleigh_drag_work_w, 0.0);
}

MPS_TEST_CASE("reported thermal power is the discrete total-energy derivative") {
  const auto parameters = config();
  const mps::DryHydrostaticDriver driver(parameters);
  const auto state = driver.initial_state();
  const auto derived = driver.diagnose(state);
  const auto forcing =
      mps::held_suarez_tendency(driver.grid(), coordinate(parameters), derived,
                                state.surface_pressure_pa, parameters.planet);
  constexpr mps::Real epsilon_s = 100.0;
  auto plus = state;
  auto minus = state;
  for (std::size_t index = 0; index < state.potential_temperature_mass_k_kg_m2.size();
       ++index) {
    plus.potential_temperature_mass_k_kg_m2[index] +=
        epsilon_s * forcing.potential_temperature_mass_k_kg_m2_s[index];
    minus.potential_temperature_mass_k_kg_m2[index] -=
        epsilon_s * forcing.potential_temperature_mass_k_kg_m2_s[index];
  }
  const auto plus_derived = driver.diagnose(plus);
  const auto minus_derived = driver.diagnose(minus);
  const auto plus_energy = mps::diagnose_dry_hydrostatic_budgets(
      driver.grid(), plus, plus_derived, parameters.planet);
  const auto minus_energy = mps::diagnose_dry_hydrostatic_budgets(
      driver.grid(), minus, minus_derived, parameters.planet);
  const mps::Real numerical_rate =
      (plus_energy.total_energy_j - minus_energy.total_energy_j) / (2.0 * epsilon_s);
  MPS_CHECK_NEAR(numerical_rate, forcing.diagnostics.thermal_energy_rate_w,
                 2e-9 * std::abs(forcing.diagnostics.thermal_energy_rate_w));
}

int main() { return mps::test::run_all(); }
