#pragma once

#include <cstdint>
#include <vector>

#include "myplanetsim/config/experiment_config.hpp"
#include "myplanetsim/dynamics/dry_hydrostatic_flux.hpp"
#include "myplanetsim/dynamics/dry_hydrostatic_state.hpp"
#include "myplanetsim/grid/cubed_sphere_grid.hpp"

namespace mps {

struct DryHydrostaticFaceStates {
  DryHydrostaticPrimitive left;
  DryHydrostaticPrimitive right;
};

struct DryHydrostaticReconstruction {
  std::size_t levels;
  std::vector<DryHydrostaticFaceStates> edge_levels;
  std::uint64_t limiter_activations;

  [[nodiscard]] const DryHydrostaticFaceStates& at(std::size_t edge,
                                                   std::size_t level) const;
};

[[nodiscard]] DryHydrostaticReconstruction reconstruct_dry_hydrostatic_face_states(
    const CubedSphereGrid& grid, const DryHydrostaticDerived& derived,
    ReconstructionKind reconstruction, LimiterKind limiter);

}  // namespace mps
