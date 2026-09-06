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

MPS_TEST_CASE("theta can use a smooth reconstruction without unbounding tracer") {
  constexpr std::size_t levels = 4;
  mps::DryHydrostaticState state{
      .surface_pressure_pa = {100000},
      .horizontal_momentum_mass_kg_m_s = {{1, 0, 0}, {2, 0, 0}, {4, 0, 0}, {8, 0, 0}},
      .potential_temperature_mass_k_kg_m2 = {1, 2, 4, 8},
      .tracer_mass_kg_m2 = {1, 2, 4, 8}};
  mps::DryHydrostaticDerived derived{
      .cells = 1, .levels = levels, .air_mass_kg_m2 = {1, 1, 1, 1}};
  mps::DryHydrostaticTransportTendency horizontal{
      .air_mass = {1, 1, -1, -1},
      .momentum = {{}, {}, {}, {}},
      .potential_temperature_mass = {0, 0, 0, 0},
      .tracer_mass = {0, 0, 0, 0}};
  const std::vector<mps::Real> b_half{0, .25, .5, .75, 1};

  mps::DryHydrostaticCoupling shared_limiter;
  mps::DryHydrostaticCouplingWorkspace shared_workspace;
  mps::couple_dry_hydrostatic_columns(
      state, derived, horizontal, b_half, 10, mps::VerticalTransportScheme::kLinear,
      mps::VerticalLimiterKind::kMinmod, shared_limiter, shared_workspace);

  mps::DryHydrostaticCoupling split_limiter;
  mps::DryHydrostaticCouplingWorkspace split_workspace;
  mps::couple_dry_hydrostatic_columns(
      state, derived, horizontal, b_half, 10, mps::VerticalTransportScheme::kLinear,
      mps::VerticalLimiterKind::kMinmod, mps::VerticalLimiterKind::kNone, split_limiter,
      split_workspace);

  bool theta_changed = false;
  for (std::size_t level = 0; level < levels; ++level) {
    theta_changed =
        theta_changed || split_limiter.tendency.potential_temperature_mass[level] !=
                             shared_limiter.tendency.potential_temperature_mass[level];
    MPS_CHECK_EQ(split_limiter.tendency.tracer_mass[level],
                 shared_limiter.tendency.tracer_mass[level]);
    MPS_CHECK_EQ(split_limiter.tendency.momentum[level].x,
                 shared_limiter.tendency.momentum[level].x);
  }
  MPS_CHECK(theta_changed);
}
int main() { return mps::test::run_all(); }
