#include <algorithm>
#include <cmath>
#include <sstream>

#include "myplanetsim/diagnostics/dry_hydrostatic_diagnostics.hpp"
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
  for (std::size_t cell = 0; cell < n; ++cell) {
    MPS_CHECK_NEAR(mps::norm(s.geopotential_gradient_kg_m_s2[cell]), 0, 1e-12);
    MPS_CHECK_NEAR(mps::norm(s.pressure_correction_kg_m_s2[cell]), 0, 1e-12);
    MPS_CHECK_NEAR(mps::norm(s.pressure_gradient_kg_m_s2[cell]), 0, 1e-12);
  }
}

namespace {
double sloping_surface_cancellation_error(const mps::Index resolution,
                                          const std::size_t levels) {
  const mps::CubedSphereGrid grid(resolution, 6.37122e6);
  mps::DryHydrostaticDerived derived{.cells = grid.cell_count(), .levels = levels};
  const auto volume = derived.cells * levels;
  derived.pressure_pa.resize(volume);
  derived.air_mass_kg_m2.assign(volume, 1000);
  derived.velocity_m_s.assign(volume, {});
  derived.temperature_k.resize(volume);
  derived.geopotential_m2_s2.resize(volume);
  constexpr double rd = 287;
  constexpr double temperature = 280;
  for (std::size_t cell = 0; cell < derived.cells; ++cell) {
    const auto position = grid.cells()[cell].center;
    const double shape = 0.35 * position.x - 0.2 * position.y + 0.1 * position.z;
    for (std::size_t level = 0; level < levels; ++level) {
      const auto offset = mps::dry_hydrostatic_offset(cell, level, levels);
      const double reference_pressure = 30000 + 5000 * static_cast<double>(level + 1);
      derived.pressure_pa[offset] = reference_pressure * std::exp(shape);
      derived.temperature_k[offset] = temperature;
      derived.geopotential_m2_s2[offset] = -rd * temperature * shape;
    }
  }
  const mps::PlanetParameters planet{grid.radius_m(), 0, 9.80616, rd, 1004.5, 100000};
  const auto sources = mps::dry_hydrostatic_sources(grid, derived, planet);
  double maximum = 0;
  for (std::size_t offset = 0; offset < volume; ++offset) {
    const auto reconstructed = sources.geopotential_gradient_kg_m_s2[offset] +
                               sources.pressure_correction_kg_m_s2[offset];
    MPS_CHECK_NEAR(mps::norm(reconstructed - sources.pressure_gradient_kg_m_s2[offset]),
                   0, 1e-12);
    maximum = std::max(maximum, mps::norm(sources.pressure_gradient_kg_m_s2[offset]));
  }
  return maximum;
}
}  // namespace

MPS_TEST_CASE("sloping hybrid pressure-gradient cancellation converges") {
  const auto coarse = sloping_surface_cancellation_error(6, 2);
  const auto medium = sloping_surface_cancellation_error(12, 2);
  const auto fine = sloping_surface_cancellation_error(24, 2);
  MPS_CHECK(medium < 0.65 * coarse);
  MPS_CHECK(fine < 0.65 * medium);
  MPS_CHECK_NEAR(sloping_surface_cancellation_error(12, 4), medium, 1e-8 * medium);
}

MPS_TEST_CASE("reference-state pressure gradient preserves its hydrostatic state") {
  const mps::CubedSphereGrid grid(8, 6.37122e6);
  mps::DryHydrostaticDerived derived{.cells = grid.cell_count(), .levels = 1};
  derived.pressure_pa.resize(derived.cells);
  derived.air_mass_kg_m2.assign(derived.cells, 1000);
  derived.velocity_m_s.assign(derived.cells, {});
  derived.temperature_k.assign(derived.cells, 280);
  derived.geopotential_m2_s2.resize(derived.cells);
  for (std::size_t cell = 0; cell < derived.cells; ++cell) {
    const double shape = 0.35 * grid.cells()[cell].center.x;
    derived.pressure_pa[cell] = 50000 * std::exp(shape);
    derived.geopotential_m2_s2[cell] = -287 * 280 * shape;
  }
  const mps::PlanetParameters planet{grid.radius_m(), 0, 9.80616, 287, 1004.5, 100000};
  const auto reference =
      mps::make_dry_hydrostatic_pressure_reference(grid, derived, planet);
  mps::DryHydrostaticSources sources;
  mps::DryHydrostaticSourcesWorkspace workspace;
  mps::dry_hydrostatic_sources(grid, derived, planet, reference, sources, workspace);
  for (const auto force : sources.pressure_gradient_kg_m_s2)
    MPS_CHECK_EQ(mps::norm(force), 0.0);

  derived.geopotential_m2_s2.front() += 1.0;
  mps::dry_hydrostatic_sources(grid, derived, planet, reference, sources, workspace);
  const bool any_nonzero =
      std::ranges::any_of(sources.pressure_gradient_kg_m_s2,
                          [](const mps::Vec3 force) { return mps::norm(force) > 0.0; });
  MPS_CHECK(any_nonzero);
}

MPS_TEST_CASE("terrain budgets and absolute pressure velocity are diagnostics") {
  mps::CubedSphereGrid grid(4, 2);
  mps::DryHydrostaticDerived derived{.cells = grid.cell_count(), .levels = 1};
  derived.pressure_pa.assign(derived.cells, 50000);
  derived.air_mass_kg_m2.assign(derived.cells, 1000);
  derived.velocity_m_s.assign(derived.cells, {1, 0, 0});
  derived.temperature_k.assign(derived.cells, 280);
  derived.geopotential_m2_s2.resize(derived.cells);
  std::vector<mps::Real> surface(derived.cells);
  for (std::size_t cell = 0; cell < derived.cells; ++cell) {
    surface[cell] = 10 * grid.cells()[cell].center.x;
    derived.geopotential_m2_s2[cell] = surface[cell];
  }
  const mps::PlanetParameters planet{2, 0, 10, 287, 1004, 100000};
  mps::DryHydrostaticState state;
  state.surface_pressure_pa.assign(derived.cells, 100000);
  const auto sources = mps::dry_hydrostatic_sources(grid, derived, planet);
  const auto diagnostics =
      mps::diagnose_terrain_budgets(grid, state, derived, sources, surface, planet);
  MPS_CHECK(diagnostics.maximum_surface_height_m > 0);
  MPS_CHECK(diagnostics.maximum_surface_slope > 0);
  MPS_CHECK(std::isfinite(diagnostics.axial_torque_residual_n_m));
  std::ostringstream output;
  mps::write_terrain_diagnostics(output, diagnostics);
  MPS_CHECK(output.str().find("terrain.pressure_work_w") != std::string::npos);
  MPS_CHECK_NEAR(
      mps::absolute_pressure_velocity_pa_s(.5, 2, {3, 0, 0}, {4, 0, 0}, 5, 10), 63, 0);
}
int main() { return mps::test::run_all(); }
