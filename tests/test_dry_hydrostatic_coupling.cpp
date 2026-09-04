#include "myplanetsim/dynamics/dry_hydrostatic_coupling.hpp"
#include "support/test.hpp"
MPS_TEST_CASE("layer convergence closes at impermeable boundaries") {
  mps::DryHydrostaticState s{.surface_pressure_pa = {100000},
                             .horizontal_momentum_mass_kg_m_s = {{}, {}, {}},
                             .potential_temperature_mass_k_kg_m2 = {3, 3, 3},
                             .tracer_mass_kg_m2 = {1, 1, 1}};
  mps::DryHydrostaticDerived d{.cells = 1, .levels = 3, .air_mass_kg_m2 = {1, 1, 1}};
  mps::DryHydrostaticTransportTendency h{.air_mass = {1, -2, 1},
                                         .momentum = {{}, {}, {}},
                                         .potential_temperature_mass = {3, -6, 3},
                                         .tracer_mass = {1, -2, 1}};
  auto out = mps::couple_dry_hydrostatic_columns(
      s, d, h, std::vector<mps::Real>{0, .2, .7, 1}, 10,
      mps::VerticalTransportScheme::kDonorCell, mps::VerticalLimiterKind::kNone);
  MPS_CHECK_NEAR(out.interface_mass_flux_kg_m2_s.front(), 0, 0);
  MPS_CHECK_NEAR(out.interface_mass_flux_kg_m2_s.back(), 0, 1e-14);
  MPS_CHECK_NEAR(out.maximum_continuity_residual_pa_s, 0, 1e-14);
}
int main() { return mps::test::run_all(); }
