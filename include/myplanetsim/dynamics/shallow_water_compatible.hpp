#pragma once

#include "myplanetsim/dynamics/shallow_water_rhs.hpp"
#include "myplanetsim/grid/dual_topology.hpp"

namespace mps {

// Vector-invariant C-grid tendency built on the compatible diagnostic
// operators. The shared unique-edge mass flux makes the depth update
// structurally conservative, and the planetary and relative vorticity enter
// through one potential vorticity flux instead of a separate Coriolis source.
//
// The scheme is centred: it reconstructs no face states, so
// `shallow_water.reconstruction` and `shallow_water.limiter` are unused.
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
