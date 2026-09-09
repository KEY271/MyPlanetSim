#pragma once

#include "myplanetsim/config/experiment_config.hpp"
#include "myplanetsim/grid/cubed_sphere_grid.hpp"

namespace mps {

// Explicit finite-volume diffusion stability bound shared by atmospheric cores.
[[nodiscard]] Real stable_diffusion_time_step(const CubedSphereGrid& grid,
                                              DiffusionKind kind,
                                              Real diffusion_coefficient,
                                              Real maximum_time_step_s);

}  // namespace mps
