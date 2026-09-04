#pragma once
#include "myplanetsim/dynamics/dry_hydrostatic_driver.hpp"
namespace mps {
[[nodiscard]] DryHydrostaticState initialize_dry_hydrostatic_benchmark(
    const ExperimentConfig& config, const CubedSphereGrid& grid,
    const AtmosphericHybridCoordinate& coordinate);
}
