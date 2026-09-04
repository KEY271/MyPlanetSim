#include "myplanetsim/dynamics/dry_hydrostatic_state.hpp"
#include "myplanetsim/grid/cubed_sphere_grid.hpp"
#include "support/test.hpp"

MPS_TEST_CASE("dry hydrostatic state is cell-major and round trips") {
  mps::DryHydrostaticState s{.surface_pressure_pa = {90000, 100000},
                             .horizontal_momentum_mass_kg_m_s = {{1, 2, 3},
                                                                 {4, 5, 6},
                                                                 {7, 8, 9},
                                                                 {10, 11, 12}},
                             .potential_temperature_mass_k_kg_m2 = {13, 14, 15, 16},
                             .tracer_mass_kg_m2 = {17, 18, 19, 20}};
  auto flat = mps::flatten_dry_hydrostatic_state(s, 2);
  auto copy = mps::unflatten_dry_hydrostatic_state(0, 0, flat, 2, 2);
  MPS_CHECK_EQ(flat.size(), 22U);
  MPS_CHECK_EQ(copy.horizontal_momentum_mass_kg_m_s[2].y, 8.0);
  MPS_CHECK_EQ(mps::dry_hydrostatic_offset(1, 0, 2), 2U);
  MPS_CHECK_THROWS_AS(mps::unflatten_dry_hydrostatic_state(0, 0, flat, 3, 2),
                      std::invalid_argument);
}

MPS_TEST_CASE("surface checkpoint layout appends one optional cell tail") {
  mps::DryHydrostaticState state{
      .time_s = 4,
      .step = 2,
      .surface_pressure_pa = {90000, 100000},
      .horizontal_momentum_mass_kg_m_s = {{1, 2, 3}, {4, 5, 6}},
      .potential_temperature_mass_k_kg_m2 = {7, 8},
      .tracer_mass_kg_m2 = {9, 10},
      .surface_temperature_k = {280, 300}};
  const auto atmospheric = mps::flatten_dry_hydrostatic_state(state, 1);
  const auto surface = mps::flatten_dry_hydrostatic_surface_state(state, 1);
  MPS_CHECK_EQ(atmospheric.size(), 12U);
  MPS_CHECK_EQ(surface.size(), 14U);
  for (std::size_t index = 0; index < atmospheric.size(); ++index)
    MPS_CHECK_EQ(surface[index], atmospheric[index]);
  const auto restored = mps::unflatten_dry_hydrostatic_surface_state(
      state.time_s, state.step, surface, 2, 1);
  MPS_CHECK(restored.surface_temperature_k == state.surface_temperature_k);
  MPS_CHECK_THROWS_AS(
      mps::unflatten_dry_hydrostatic_surface_state(0, 0, atmospheric, 2, 1),
      std::invalid_argument);
}

MPS_TEST_CASE("cellwise lower boundary shifts only diagnosed geopotential") {
  const mps::AtmosphericHybridCoordinate coordinate({{1000, 0}, {0, 1}}, 90000, 110000,
                                                    100);
  const mps::PlanetParameters planet{2, 0, 10, 287, 1004, 100000};
  mps::DryHydrostaticState state{
      .surface_pressure_pa = {100000, 100000},
      .horizontal_momentum_mass_kg_m_s = {{1, 2, 0}, {1, 2, 0}},
      .potential_temperature_mass_k_kg_m2 = {3000000, 3000000},
      .tracer_mass_kg_m2 = {100, 100}};
  const auto flat = mps::diagnose_dry_hydrostatic_state(state, coordinate, planet);
  const std::vector<mps::Real> lower_boundary{123, 456};
  const auto terrain =
      mps::diagnose_dry_hydrostatic_state(state, coordinate, planet, lower_boundary);
  for (std::size_t cell = 0; cell < 2; ++cell) {
    MPS_CHECK_EQ(terrain.pressure_pa[cell], flat.pressure_pa[cell]);
    MPS_CHECK_EQ(terrain.temperature_k[cell], flat.temperature_k[cell]);
    MPS_CHECK_EQ(terrain.velocity_m_s[cell].x, flat.velocity_m_s[cell].x);
    MPS_CHECK_NEAR(terrain.geopotential_m2_s2[cell] - flat.geopotential_m2_s2[cell],
                   lower_boundary[cell], 1e-12);
  }
  MPS_CHECK_THROWS_AS(mps::diagnose_dry_hydrostatic_state(state, coordinate, planet,
                                                          std::vector<mps::Real>{1}),
                      std::invalid_argument);
}
int main() { return mps::test::run_all(); }
