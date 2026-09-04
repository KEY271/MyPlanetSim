#include "myplanetsim/dynamics/dry_hydrostatic_sources.hpp"

#include <stdexcept>

#include "myplanetsim/numerics/spherical_operators.hpp"
namespace mps {
DryHydrostaticSources dry_hydrostatic_sources(const CubedSphereGrid& grid,
                                              const DryHydrostaticDerived& d,
                                              const PlanetParameters& p) {
  if (d.cells != grid.cell_count() || d.levels == 0)
    throw std::invalid_argument("dry source shape mismatch");
  DryHydrostaticSources out;
  out.pressure_gradient_kg_m_s2.resize(d.cells * d.levels);
  out.coriolis_kg_m_s2.resize(d.cells * d.levels);
  const Vec3 omega{0, 0, p.rotation_rate_rad_s};
  std::vector<Real> pressure(d.cells), phi(d.cells);
  for (std::size_t k = 0; k < d.levels; ++k) {
    for (std::size_t c = 0; c < d.cells; ++c) {
      auto n = dry_hydrostatic_offset(c, k, d.levels);
      pressure[c] = d.pressure_pa[n];
      phi[c] = d.geopotential_m2_s2[n];
    }
    auto gp = least_squares_gradient(grid, pressure);
    auto gf = least_squares_gradient(grid, phi);
    for (std::size_t c = 0; c < d.cells; ++c) {
      auto n = dry_hydrostatic_offset(c, k, d.levels);
      auto centre = grid.cells()[c].center;
      auto alpha = p.gas_constant_j_kg_k * d.temperature_k[n] / d.pressure_pa[n];
      out.pressure_gradient_kg_m_s2[n] =
          -d.air_mass_kg_m2[n] * project_tangent(gf[c] + alpha * gp[c], centre);
      out.coriolis_kg_m_s2[n] =
          -2 * d.air_mass_kg_m2[n] *
          project_tangent(cross(omega, d.velocity_m_s[n]), centre);
    }
  }
  return out;
}
}  // namespace mps
