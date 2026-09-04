#include <cmath>
#include <stdexcept>
#include <vector>

#include "myplanetsim/dynamics/dry_hydrostatic_diffusion.hpp"
#include "myplanetsim/dynamics/dry_hydrostatic_driver.hpp"
#include "myplanetsim/numerics/spherical_operators.hpp"
#include "support/test.hpp"

// Gates for the dry horizontal diffusion that ADR 0006 specified and nothing
// implemented until Phase 9, leaving the ADR 0009 six-run matrix blocked. Before the
// wiring the configuration parsed and was then discarded, so "diffusion changes the
// run" and "diffusion is reported separately" both failed while "none is unchanged"
// passed vacuously.

namespace {

mps::ExperimentConfig config(const mps::DiffusionKind kind,
                             const mps::Real coefficient) {
  mps::ExperimentConfig result{
      .kind = mps::ExperimentKind::kDryHydrostatic,
      .planet = {6371220, 7.29212e-5, 9.80616, 287, 1004, 100000},
      .run = {0, 600, 300, 0},
      .grid = {4},
      .vertical = {.levels = 4,
                   .a_half_pa = {1000, 750, 500, 250, 0},
                   .b_half = {0, .25, .5, .75, 1},
                   .surface_pressure_pa = 100000,
                   .minimum_surface_pressure_pa = 50000,
                   .maximum_surface_pressure_pa = 150000,
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
  result.dry_hydrostatic.test_case = mps::DryHydrostaticTestCase::kHeldSuarez;
  result.dry_hydrostatic.diffusion_kind = kind;
  result.dry_hydrostatic.diffusion_coefficient = coefficient;
  return result;
}

[[nodiscard]] mps::DryHydrostaticState advanced(const mps::DiffusionKind kind,
                                                const mps::Real coefficient,
                                                mps::Real& diffusion_energy_j) {
  const mps::DryHydrostaticDriver driver(config(kind, coefficient));
  auto state = driver.initial_state();
  diffusion_energy_j = 0.0;
  driver.advance(state, 600.0,
                 [&](const mps::DryHydrostaticState&, const mps::DryHydrostaticDerived*,
                     const mps::DryHydrostaticStepDiagnostics& step) {
                   diffusion_energy_j += step.diffusion_energy_contribution_j;
                 });
  return state;
}

[[nodiscard]] mps::Real enstrophy(const mps::DryHydrostaticDriver& driver,
                                  const mps::DryHydrostaticState& state) {
  const auto derived = driver.diagnose(state);
  const auto& grid = driver.grid();
  mps::Real total = 0.0;
  std::vector<mps::Vec3> level(derived.cells);
  for (std::size_t k = 0; k < derived.levels; ++k) {
    for (std::size_t c = 0; c < derived.cells; ++c) {
      level[c] =
          derived.velocity_m_s[mps::dry_hydrostatic_offset(c, k, derived.levels)];
    }
    const auto vorticity = mps::finite_volume_curl(grid, level);
    for (std::size_t c = 0; c < derived.cells; ++c) {
      total += grid.cells()[c].area_m2 * vorticity[c] * vorticity[c];
    }
  }
  return total;
}

}  // namespace

MPS_TEST_CASE("disabled dry diffusion contributes exactly zero") {
  // Every shipped preset selects `none`, so this exact zero is what keeps the Phase
  // 5--8 regression baselines byte-exact across the wiring.
  // `phase5.dry_hydrostatic_baseline` checks the resulting fingerprints; this checks
  // the reason they cannot move.
  const mps::DryHydrostaticDriver driver(config(mps::DiffusionKind::kNone, 0.0));
  const auto derived = driver.diagnose(driver.initial_state());
  const auto tendency = mps::dry_hydrostatic_diffusion_tendency(
      driver.grid(), derived, mps::DiffusionKind::kNone, 0.0);
  MPS_CHECK_EQ(tendency.kinetic_energy_rate_w, 0.0);
  for (std::size_t n = 0; n < tendency.potential_temperature_mass.size(); ++n) {
    MPS_CHECK_EQ(mps::norm(tendency.momentum[n]), 0.0);
    MPS_CHECK_EQ(tendency.potential_temperature_mass[n], 0.0);
    MPS_CHECK_EQ(tendency.tracer_mass[n], 0.0);
  }

  mps::Real reported = 0.0;
  static_cast<void>(advanced(mps::DiffusionKind::kNone, 0.0, reported));
  MPS_CHECK_EQ(reported, 0.0);
}

MPS_TEST_CASE("dry diffusion damps vorticity without moving mass") {
  const mps::DryHydrostaticDriver plain(config(mps::DiffusionKind::kNone, 0.0));
  const mps::DryHydrostaticDriver diffusive(
      config(mps::DiffusionKind::kLaplacian, 1.0e5));
  mps::Real ignored = 0.0;
  mps::Real dissipated = 0.0;
  const auto undiffused = advanced(mps::DiffusionKind::kNone, 0.0, ignored);
  const auto diffused = advanced(mps::DiffusionKind::kLaplacian, 1.0e5, dissipated);

  // Diffusion is a dissipative term: it removes enstrophy and kinetic energy.
  MPS_CHECK(enstrophy(diffusive, diffused) < enstrophy(plain, undiffused));
  MPS_CHECK(dissipated < 0.0);

  // Mass and surface pressure are untouched, so the two runs still agree on them.
  const auto plain_derived = plain.diagnose(undiffused);
  const auto diffusive_derived = diffusive.diagnose(diffused);
  for (std::size_t c = 0; c < undiffused.surface_pressure_pa.size(); ++c) {
    MPS_CHECK_NEAR(diffused.surface_pressure_pa[c], undiffused.surface_pressure_pa[c],
                   1.0e-9 * undiffused.surface_pressure_pa[c]);
  }
  for (std::size_t n = 0; n < plain_derived.air_mass_kg_m2.size(); ++n) {
    MPS_CHECK_NEAR(diffusive_derived.air_mass_kg_m2[n], plain_derived.air_mass_kg_m2[n],
                   1.0e-9 * plain_derived.air_mass_kg_m2[n]);
  }
}

MPS_TEST_CASE("an inactive kind and a nonzero coefficient are mutually exclusive") {
  const mps::DryHydrostaticDriver driver(config(mps::DiffusionKind::kNone, 0.0));
  const auto derived = driver.diagnose(driver.initial_state());
  MPS_CHECK_THROWS_AS(mps::dry_hydrostatic_diffusion_tendency(
                          driver.grid(), derived, mps::DiffusionKind::kNone, 1.0),
                      std::invalid_argument);
  MPS_CHECK_THROWS_AS(mps::dry_hydrostatic_diffusion_tendency(
                          driver.grid(), derived, mps::DiffusionKind::kLaplacian, 0.0),
                      std::invalid_argument);
}

int main() { return mps::test::run_all(); }
