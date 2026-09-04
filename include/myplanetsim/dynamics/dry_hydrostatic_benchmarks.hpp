#pragma once
#include "myplanetsim/dynamics/dry_hydrostatic_state.hpp"
#include "myplanetsim/dynamics/surface_orography.hpp"
namespace mps {
[[nodiscard]] DryHydrostaticState initialize_dry_hydrostatic_benchmark(
    const ExperimentConfig& config, const CubedSphereGrid& grid,
    const AtmosphericHybridCoordinate& coordinate,
    const SurfaceOrography& orography);
}
