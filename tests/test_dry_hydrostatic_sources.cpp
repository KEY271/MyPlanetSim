#include "myplanetsim/dynamics/dry_hydrostatic_sources.hpp"
#include "support/test.hpp"
MPS_TEST_CASE("uniform hydrostatic fields have zero pressure source") {
  mps::CubedSphereGrid g(2, 2);
  mps::DryHydrostaticDerived d{.cells = g.cell_count(), .levels = 1};
  auto n = g.cell_count();
  d.pressure_pa.assign(n, 50000);
  d.air_mass_kg_m2.assign(n, 1000);
  d.velocity_m_s.assign(n, {});
  d.temperature_k.assign(n, 280);
  d.geopotential_m2_s2.assign(n, 10);
  mps::PlanetParameters p{2, 0.1, 10, 287, 1004, 100000};
  auto s = mps::dry_hydrostatic_sources(g, d, p);
  for (auto v : s.pressure_gradient_kg_m_s2) MPS_CHECK_NEAR(mps::norm(v), 0, 1e-12);
}
int main() { return mps::test::run_all(); }
