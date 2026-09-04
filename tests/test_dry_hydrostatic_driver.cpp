#include "myplanetsim/dynamics/dry_hydrostatic_driver.hpp"
#include "support/test.hpp"
MPS_TEST_CASE("uniform isothermal rest remains stationary") {
  mps::ExperimentConfig c{
      .kind = mps::ExperimentKind::kDryHydrostatic,
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
int main() { return mps::test::run_all(); }
