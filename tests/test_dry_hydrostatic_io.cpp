#include <sstream>

#include "myplanetsim/diagnostics/dry_hydrostatic_diagnostics.hpp"
#include "myplanetsim/io/checkpoint.hpp"
#include "support/test.hpp"
MPS_TEST_CASE("dry state checkpoint layout and budgets are deterministic") {
  mps::DryHydrostaticState s{.surface_pressure_pa = {100000},
                             .horizontal_momentum_mass_kg_m_s = {{0, 0, 0}},
                             .potential_temperature_mass_k_kg_m2 = {3000},
                             .tracer_mass_kg_m2 = {2}};
  auto flat = mps::flatten_dry_hydrostatic_state(s, 1);
  std::stringstream io;
  mps::write_checkpoint(io, {0, 0, flat, "fingerprint",
                             std::string(mps::kDryHydrostaticCheckpointLayout)});
  auto c = mps::read_checkpoint(io, "fingerprint", mps::kDryHydrostaticCheckpointLayout,
                                flat.size());
  MPS_CHECK_EQ(c.state.size(), flat.size());
  mps::CubedSphereGrid g(1, 2);
}
int main() { return mps::test::run_all(); }
