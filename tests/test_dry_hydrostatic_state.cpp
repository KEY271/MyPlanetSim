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
int main() { return mps::test::run_all(); }
