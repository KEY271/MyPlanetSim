#pragma once

#include <span>
#include <vector>

#include "myplanetsim/grid/dual_topology.hpp"

namespace mps {

// Diagnostic operators of the compatible (C-grid) candidate discretization.
// The primal cell carries the depth, the primal edge carries the normal
// velocity, and the dual cell around a primal vertex carries the vorticity.
//
// The cubed-sphere dual edge is not parallel to the primal edge normal, so the
// dual-cell circulation needs both velocity components at the edge. The normal
// component is prognostic and the tangential one is reconstructed from the
// surrounding normal components.

struct CompatibleEdgeMetric {
  // Primal edge length divided by dual edge length: the diagonal Hodge star.
  Real hodge_ratio;
  // Dual edge direction projected on the edge normal and on the edge tangent.
  Real normal_projection;
  Real tangent_projection;
};

[[nodiscard]] std::vector<CompatibleEdgeMetric> compatible_edge_metrics(
    const CubedSphereGrid& grid, const CubedSphereDualTopology& dual);
[[nodiscard]] std::vector<Real> edge_normal_velocity(
    const CubedSphereGrid& grid, const CubedSphereDualTopology& dual,
    std::span<const Vec3> cell_velocity);
// Perot reconstruction of the cell tangent velocity from edge normal velocities.
[[nodiscard]] std::vector<Vec3> perot_cell_velocity(
    const CubedSphereGrid& grid, std::span<const Real> edge_normal_velocity);
[[nodiscard]] std::vector<Real> vertex_relative_vorticity(
    const CubedSphereGrid& grid, const CubedSphereDualTopology& dual,
    std::span<const Real> edge_normal_velocity);
[[nodiscard]] std::vector<Real> vertex_depth(const CubedSphereGrid& grid,
                                             const CubedSphereDualTopology& dual,
                                             std::span<const Real> cell_depth);
[[nodiscard]] std::vector<Real> vertex_potential_vorticity(
    const CubedSphereGrid& grid, const CubedSphereDualTopology& dual,
    std::span<const Real> cell_depth, std::span<const Real> edge_normal_velocity,
    Vec3 rotation_vector_rad_s);
[[nodiscard]] std::vector<Real> edge_potential_vorticity(
    const CubedSphereGrid& grid, std::span<const Real> vertex_potential_vorticity);
[[nodiscard]] std::vector<Real> cell_kinetic_energy(
    const CubedSphereGrid& grid, std::span<const Real> edge_normal_velocity);

}  // namespace mps
