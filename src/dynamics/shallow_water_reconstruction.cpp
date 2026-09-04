#include "myplanetsim/dynamics/shallow_water_reconstruction.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace mps {
namespace {

[[nodiscard]] Real barth_factor(const Real center, const Real increment,
                                const Real minimum, const Real maximum) {
  if (increment > 0.0) {
    return std::min(1.0, (maximum - center) / increment);
  }
  if (increment < 0.0) {
    return std::min(1.0, (minimum - center) / increment);
  }
  return 1.0;
}

[[nodiscard]] const CachedCellEdgeGeometry& cached_cell_edge(
    const CubedSphereGrid& grid, const std::size_t cell, const std::size_t edge) {
  for (const auto& cached : grid.cell_cache()[cell].edges) {
    if (cached.edge == edge) return cached;
  }
  throw std::logic_error("cached edge is not incident to cell");
}

[[nodiscard]] Real speed_factor(const Vec3 center, const Vec3 increment,
                                const Real maximum_speed) {
  if (norm(center + increment) <= maximum_speed) {
    return 1.0;
  }
  Real low = 0.0;
  Real high = 1.0;
  for (int iteration = 0; iteration < 54; ++iteration) {
    const Real middle = 0.5 * (low + high);
    if (norm(center + middle * increment) <= maximum_speed) {
      low = middle;
    } else {
      high = middle;
    }
  }
  return low;
}

}  // namespace

ShallowWaterReconstruction reconstruct_shallow_water_face_states(
    const CubedSphereGrid& grid, const ShallowWaterState& state,
    const ReconstructionKind reconstruction, const LimiterKind limiter,
    const Real depth_floor_m) {
  validate_shallow_water_state(grid, state, depth_floor_m);
  ShallowWaterReconstruction result{
      .edges = std::vector<ShallowWaterFaceStates>(grid.edge_count()),
      .limiter_activations = 0,
  };
  std::vector<Vec3> velocity(grid.cell_count());
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
    velocity[cell] = state.velocity(cell);
  }
  if (reconstruction == ReconstructionKind::kPiecewiseConstant) {
    for (const auto& edge : grid.edges()) {
      const auto& cached = grid.edge_cache()[edge.id];
      const std::size_t left = cached.left_cell;
      const std::size_t right = cached.right_cell;
      result.edges[edge.id] = {
          .left = {.depth_m = state.depth[left],
                   .velocity_m_s = project_tangent(velocity[left], edge.center)},
          .right = {.depth_m = state.depth[right],
                    .velocity_m_s = project_tangent(velocity[right], edge.center)},
      };
    }
    return result;
  }

  const auto depth_gradient = least_squares_gradient(grid, state.depth);
  const auto velocity_gradient = least_squares_vector_gradient(grid, velocity);
  std::vector<Real> depth_factor(grid.cell_count(), 1.0);
  std::vector<Real> velocity_factor(grid.cell_count(), 1.0);
  if (limiter == LimiterKind::kBarthJespersen) {
    for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
      Real minimum_depth = state.depth[cell];
      Real maximum_depth = state.depth[cell];
      Real maximum_speed = norm(velocity[cell]);
      for (const auto& edge : grid.cell_cache()[cell].edges) {
        minimum_depth = std::min(minimum_depth, state.depth[edge.neighbor]);
        maximum_depth = std::max(maximum_depth, state.depth[edge.neighbor]);
        maximum_speed = std::max(maximum_speed, norm(velocity[edge.neighbor]));
      }
      const Real positive_floor =
          std::nextafter(depth_floor_m, std::numeric_limits<Real>::infinity());
      for (const auto& cached_edge : grid.cell_cache()[cell].edges) {
        const auto& edge = grid.edges()[cached_edge.edge];
        const Real depth_increment =
            dot(depth_gradient[cell], cached_edge.face_displacement_m);
        depth_factor[cell] = std::min(depth_factor[cell],
                                      barth_factor(state.depth[cell], depth_increment,
                                                   minimum_depth, maximum_depth));
        if (depth_increment < 0.0) {
          depth_factor[cell] =
              std::min(depth_factor[cell],
                       (state.depth[cell] - positive_floor) / -depth_increment);
        }

        const auto& edge_cache = grid.edge_cache()[cached_edge.edge];
        const EdgeTangentBasis basis{edge_cache.normal, edge_cache.tangent};
        const Vec3 center_at_face = project_tangent(velocity[cell], edge.center);
        const Vec3 unlimited = reconstruct_tangent_vector_cached(
            grid, cell, velocity[cell], edge.center,
            cached_edge.normalized_face_displacement_m, velocity_gradient[cell]);
        const Vec3 increment = unlimited - center_at_face;
        Real minimum_normal = dot(center_at_face, basis.normal);
        Real maximum_normal = minimum_normal;
        Real minimum_tangent = dot(center_at_face, basis.tangent);
        Real maximum_tangent = minimum_tangent;
        for (const auto& neighbor_edge : grid.cell_cache()[cell].edges) {
          const Vec3 neighbor_at_face =
              project_tangent(velocity[neighbor_edge.neighbor], edge.center);
          const Real normal = dot(neighbor_at_face, basis.normal);
          const Real tangent = dot(neighbor_at_face, basis.tangent);
          minimum_normal = std::min(minimum_normal, normal);
          maximum_normal = std::max(maximum_normal, normal);
          minimum_tangent = std::min(minimum_tangent, tangent);
          maximum_tangent = std::max(maximum_tangent, tangent);
        }
        velocity_factor[cell] = std::min(
            velocity_factor[cell],
            barth_factor(dot(center_at_face, basis.normal),
                         dot(increment, basis.normal), minimum_normal, maximum_normal));
        velocity_factor[cell] = std::min(
            velocity_factor[cell], barth_factor(dot(center_at_face, basis.tangent),
                                                dot(increment, basis.tangent),
                                                minimum_tangent, maximum_tangent));
        velocity_factor[cell] =
            std::min(velocity_factor[cell],
                     speed_factor(center_at_face, increment, maximum_speed));
      }
      depth_factor[cell] = std::clamp(depth_factor[cell], 0.0, 1.0);
      velocity_factor[cell] = std::clamp(velocity_factor[cell], 0.0, 1.0);
      if (depth_factor[cell] < 1.0 - 1.0e-14 || velocity_factor[cell] < 1.0 - 1.0e-14) {
        ++result.limiter_activations;
      }
    }
  }

  for (const auto& edge : grid.edges()) {
    const auto& edge_cache = grid.edge_cache()[edge.id];
    const std::size_t left = edge_cache.left_cell;
    const std::size_t right = edge_cache.right_cell;
    const auto reconstruct = [&](const std::size_t cell) {
      const auto& cell_edge = cached_cell_edge(grid, cell, edge.id);
      const Real depth =
          state.depth[cell] +
          depth_factor[cell] * dot(depth_gradient[cell], cell_edge.face_displacement_m);
      const Vec3 center_velocity = project_tangent(velocity[cell], edge.center);
      const Vec3 unlimited = reconstruct_tangent_vector_cached(
          grid, cell, velocity[cell], edge.center,
          cell_edge.normalized_face_displacement_m, velocity_gradient[cell]);
      return ShallowWaterPrimitive{
          .depth_m = depth,
          .velocity_m_s =
              center_velocity + velocity_factor[cell] * (unlimited - center_velocity),
      };
    };
    result.edges[edge.id] = {.left = reconstruct(left), .right = reconstruct(right)};
  }
  return result;
}

}  // namespace mps
