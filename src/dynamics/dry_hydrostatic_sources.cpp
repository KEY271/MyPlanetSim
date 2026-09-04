#include "myplanetsim/dynamics/dry_hydrostatic_sources.hpp"

#include <stdexcept>

#include "myplanetsim/numerics/spherical_operators.hpp"
namespace mps {
DryHydrostaticSources dry_hydrostatic_sources(const CubedSphereGrid& grid,
                                              const DryHydrostaticDerived& d,
                                              const PlanetParameters& p) {
  DryHydrostaticSources result;
  DryHydrostaticSourcesWorkspace workspace;
  dry_hydrostatic_sources(grid, d, p, result, workspace);
  return result;
}

void dry_hydrostatic_sources(const CubedSphereGrid& grid,
                             const DryHydrostaticDerived& d, const PlanetParameters& p,
                             DryHydrostaticSources& out,
                             DryHydrostaticSourcesWorkspace& workspace) {
  if (d.cells != grid.cell_count() || d.levels == 0)
    throw std::invalid_argument("dry source shape mismatch");
  out.geopotential_gradient_kg_m_s2.resize(d.cells * d.levels);
  out.pressure_correction_kg_m_s2.resize(d.cells * d.levels);
  out.pressure_gradient_kg_m_s2.resize(d.cells * d.levels);
  out.coriolis_kg_m_s2.resize(d.cells * d.levels);
  const Vec3 omega{0, 0, p.rotation_rate_rad_s};
  workspace.pressure.resize(d.cells);
  workspace.geopotential.resize(d.cells);
  workspace.pressure_gradient.resize(d.cells);
  workspace.geopotential_gradient.resize(d.cells);
  for (std::size_t k = 0; k < d.levels; ++k) {
    for (std::size_t c = 0; c < d.cells; ++c) {
      auto n = dry_hydrostatic_offset(c, k, d.levels);
      workspace.pressure[c] = d.pressure_pa[n];
      workspace.geopotential[c] = d.geopotential_m2_s2[n];
    }
    least_squares_gradient(grid, workspace.pressure, workspace.pressure_gradient);
    least_squares_gradient(grid, workspace.geopotential,
                           workspace.geopotential_gradient);
    for (std::size_t c = 0; c < d.cells; ++c) {
      auto n = dry_hydrostatic_offset(c, k, d.levels);
      auto centre = grid.cells()[c].center;
      auto alpha = p.gas_constant_j_kg_k * d.temperature_k[n] / d.pressure_pa[n];
      out.geopotential_gradient_kg_m_s2[n] =
          -d.air_mass_kg_m2[n] *
          project_tangent(workspace.geopotential_gradient[c], centre);
      out.pressure_correction_kg_m_s2[n] =
          -d.air_mass_kg_m2[n] *
          project_tangent(alpha * workspace.pressure_gradient[c], centre);
      out.pressure_gradient_kg_m_s2[n] =
          out.geopotential_gradient_kg_m_s2[n] + out.pressure_correction_kg_m_s2[n];
      out.coriolis_kg_m_s2[n] =
          -2 * d.air_mass_kg_m2[n] *
          project_tangent(cross(omega, d.velocity_m_s[n]), centre);
    }
  }
}
}  // namespace mps
