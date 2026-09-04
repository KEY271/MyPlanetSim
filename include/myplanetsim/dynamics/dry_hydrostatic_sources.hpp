#pragma once
#include "myplanetsim/dynamics/dry_hydrostatic_state.hpp"
#include "myplanetsim/grid/cubed_sphere_grid.hpp"
namespace mps {
struct DryHydrostaticSources { std::vector<Vec3> pressure_gradient_kg_m_s2; std::vector<Vec3> coriolis_kg_m_s2; };
[[nodiscard]] DryHydrostaticSources dry_hydrostatic_sources(const CubedSphereGrid& grid,const DryHydrostaticDerived& derived,const PlanetParameters& planet);
}
