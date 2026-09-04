#pragma once
#include "myplanetsim/dynamics/dry_hydrostatic_state.hpp"
#include "myplanetsim/grid/cubed_sphere_grid.hpp"
namespace mps {
struct DryHydrostaticSources {
  std::vector<Vec3> geopotential_gradient_kg_m_s2;
  std::vector<Vec3> pressure_correction_kg_m_s2;
  std::vector<Vec3> pressure_gradient_kg_m_s2;
  std::vector<Vec3> coriolis_kg_m_s2;
};
struct DryHydrostaticSourcesWorkspace {
  std::vector<Real> pressure;
  std::vector<Real> geopotential;
  std::vector<Vec3> pressure_gradient;
  std::vector<Vec3> geopotential_gradient;
};
[[nodiscard]] DryHydrostaticSources dry_hydrostatic_sources(
    const CubedSphereGrid& grid, const DryHydrostaticDerived& derived,
    const PlanetParameters& planet);
void dry_hydrostatic_sources(const CubedSphereGrid& grid,
                             const DryHydrostaticDerived& derived,
                             const PlanetParameters& planet,
                             DryHydrostaticSources& result,
                             DryHydrostaticSourcesWorkspace& workspace);
}  // namespace mps
