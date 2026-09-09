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
  std::size_t tracer_count = 1;
  std::vector<DryHydrostaticFaceStates> edge_levels;
  std::vector<Real> left_tracer_mixing_ratio;
  std::vector<Real> right_tracer_mixing_ratio;
  std::uint64_t limiter_activations;

  [[nodiscard]] const DryHydrostaticFaceStates& at(std::size_t edge,
                                                   std::size_t level) const;
  [[nodiscard]] Real tracer_at(bool left, std::size_t tracer, std::size_t edge,
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

struct DryHydrostaticPreparedScalarReconstruction {
  std::vector<Vec3> gradient;
  std::vector<Real> factor;
};

struct DryHydrostaticPreparedReconstruction {
  std::size_t cells = 0;
  std::size_t levels = 0;
  std::size_t tracer_count = 0;
  ReconstructionKind kind = ReconstructionKind::kPiecewiseConstant;
  bool reconstruct_temperature = true;
  DryHydrostaticPreparedScalarReconstruction air_mass;
  DryHydrostaticPreparedScalarReconstruction potential_temperature;
  DryHydrostaticPreparedScalarReconstruction temperature;
  DryHydrostaticPreparedScalarReconstruction tracer;
  std::vector<TangentVectorGradient> velocity_gradient;
  std::vector<Real> velocity_factor;
  std::vector<std::uint8_t> tracer_is_constant;
  std::uint64_t limiter_activations = 0;
};

void prepare_dry_hydrostatic_reconstruction(
    const CubedSphereGrid& grid, const DryHydrostaticDerived& derived,
    ReconstructionKind reconstruction, LimiterKind limiter,
    DryHydrostaticPreparedReconstruction& result,
    DryHydrostaticReconstructionWorkspace& workspace,
    bool reconstruct_temperature = true);
[[nodiscard]] DryHydrostaticFaceStates reconstruct_dry_hydrostatic_edge(
    const CubedSphereGrid& grid, const DryHydrostaticDerived& derived,
    const DryHydrostaticPreparedReconstruction& prepared, std::size_t edge,
    std::size_t level);
[[nodiscard]] Real reconstruct_dry_hydrostatic_tracer_face(
    const CubedSphereGrid& grid, const DryHydrostaticDerived& derived,
    const DryHydrostaticPreparedReconstruction& prepared, bool left, std::size_t tracer,
    std::size_t edge, std::size_t level);

[[nodiscard]] DryHydrostaticReconstruction reconstruct_dry_hydrostatic_face_states(
    const CubedSphereGrid& grid, const DryHydrostaticDerived& derived,
    ReconstructionKind reconstruction, LimiterKind limiter);
void reconstruct_dry_hydrostatic_face_states(
    const CubedSphereGrid& grid, const DryHydrostaticDerived& derived,
    ReconstructionKind reconstruction, LimiterKind limiter,
    DryHydrostaticReconstruction& result,
    DryHydrostaticReconstructionWorkspace& workspace,
    bool reconstruct_temperature = true);

}  // namespace mps
