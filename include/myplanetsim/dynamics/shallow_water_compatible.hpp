#pragma once

#include "myplanetsim/dynamics/shallow_water_rhs.hpp"
#include "myplanetsim/grid/dual_topology.hpp"

namespace mps {

// Vector-invariant C-grid tendency built on the compatible diagnostic
// operators. The shared unique-edge mass flux makes the depth update
// structurally conservative; the vorticity term is written as a potential
// vorticity flux so that a constant potential vorticity stays constant.
//
// The reported components follow the reference scheme where they can:
// `flux` holds the mass divergence together with the kinetic-energy gradient,
// `pressure` the free-surface gradient, and `coriolis` the full potential
// vorticity flux, which unlike the reference scheme also carries the relative
// vorticity.
[[nodiscard]] ShallowWaterRhsComponents assemble_compatible_shallow_water_rhs(
    const CubedSphereGrid& grid, const CubedSphereDualTopology& dual,
    const ShallowWaterState& state, const ShallowWaterParameters& parameters,
    Real gravity_m_s2, Vec3 rotation_vector_rad_s);

}  // namespace mps
