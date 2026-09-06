#pragma once

#include <array>
#include <cstdint>
#include <vector>

#include "myplanetsim/config/experiment_config.hpp"
#include "myplanetsim/dynamics/dry_hydrostatic_flux.hpp"
#include "myplanetsim/dynamics/dry_hydrostatic_state.hpp"
#include "myplanetsim/grid/cubed_sphere_grid.hpp"
#include "myplanetsim/numerics/spherical_operators.hpp"

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

struct DryHydrostaticReconstructionLevelWorkspace {
  std::array<std::vector<Vec3>, 4> scalar_gradients;
  std::array<std::vector<Real>, 5> limiter_factors;
  std::vector<TangentVectorGradient> velocity_gradient;
  std::vector<Real> mass;
  std::vector<Real> potential_temperature;
  std::vector<Real> tracer;
  std::vector<Real> temperature;
  std::vector<Vec3> velocity;
};

struct DryHydrostaticReconstructionWorkspace {
  std::vector<DryHydrostaticReconstructionLevelWorkspace> workers;
};

[[nodiscard]] DryHydrostaticReconstruction reconstruct_dry_hydrostatic_face_states(
    const CubedSphereGrid& grid, const DryHydrostaticDerived& derived,
    ReconstructionKind reconstruction, LimiterKind limiter);
void reconstruct_dry_hydrostatic_face_states(
    const CubedSphereGrid& grid, const DryHydrostaticDerived& derived,
    ReconstructionKind reconstruction, LimiterKind limiter,
    DryHydrostaticReconstruction& result,
    DryHydrostaticReconstructionWorkspace& workspace);

}  // namespace mps
