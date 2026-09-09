#pragma once
#include <exception>

#include "myplanetsim/dynamics/dry_hydrostatic_state.hpp"
#include "myplanetsim/grid/cubed_sphere_grid.hpp"
namespace mps {
struct DryHydrostaticPressureReference {
  std::size_t cells = 0;
  std::size_t levels = 0;
  std::vector<Real> pressure_pa;
  std::vector<Real> geopotential_m2_s2;
  std::vector<Real> specific_volume_m3_kg;
  std::vector<Vec3> pressure_gradient_pa_m;
};
struct DryHydrostaticSources {
  std::vector<Vec3> geopotential_gradient_kg_m_s2;
  std::vector<Vec3> pressure_correction_kg_m_s2;
  std::vector<Vec3> pressure_gradient_kg_m_s2;
  std::vector<Vec3> coriolis_kg_m_s2;
};
struct DryHydrostaticSourcesLevelWorkspace {
  std::vector<Real> pressure;
  std::vector<Real> geopotential;
  std::vector<Vec3> pressure_gradient;
  std::vector<Vec3> geopotential_gradient;
};
struct DryHydrostaticSourcesWorkspace {
  std::vector<DryHydrostaticSourcesLevelWorkspace> workers;
  std::vector<std::exception_ptr> failures;
};
[[nodiscard]] DryHydrostaticSources dry_hydrostatic_sources(
    const CubedSphereGrid& grid, const DryHydrostaticDerived& derived,
    const PlanetParameters& planet);
[[nodiscard]] DryHydrostaticPressureReference make_dry_hydrostatic_pressure_reference(
    const CubedSphereGrid& grid, const DryHydrostaticDerived& derived,
    const PlanetParameters& planet);
void dry_hydrostatic_sources(const CubedSphereGrid& grid,
                             const DryHydrostaticDerived& derived,
                             const PlanetParameters& planet,
                             DryHydrostaticSources& result,
                             DryHydrostaticSourcesWorkspace& workspace);
void dry_hydrostatic_sources(const CubedSphereGrid& grid,
                             const DryHydrostaticDerived& derived,
                             const PlanetParameters& planet,
                             const DryHydrostaticPressureReference& reference,
                             DryHydrostaticSources& result,
                             DryHydrostaticSourcesWorkspace& workspace);
}  // namespace mps
