#pragma once

#include <vector>

#include "myplanetsim/config/experiment_config.hpp"
#include "myplanetsim/dynamics/dry_hydrostatic_state.hpp"
#include "myplanetsim/grid/cubed_sphere_grid.hpp"

namespace mps {

// Explicit horizontal diffusion for the dry hydrostatic core (ADR 0006, wired in
// Phase 9). The existing horizontal Laplacian and biharmonic operators are applied
// level by level to the intensive fields and then weighted by the unchanged layer mass,
// so diffusion moves momentum, heat, and tracer without moving mass. Air mass, surface
// pressure, and terrain are untouched by construction.
struct DryHydrostaticDiffusionTendency {
  std::vector<Vec3> momentum;
  std::vector<Real> potential_temperature_mass;
  std::vector<Real> tracer_mass;
  // Diffusion work on the resolved flow, negative when diffusion removes kinetic
  // energy. Reported separately so dissipation is attributed rather than hidden.
  Real kinetic_energy_rate_w = 0.0;
};

[[nodiscard]] DryHydrostaticDiffusionTendency dry_hydrostatic_diffusion_tendency(
    const CubedSphereGrid& grid, const DryHydrostaticDerived& derived,
    DiffusionKind kind, Real diffusion_coefficient);

}  // namespace mps
