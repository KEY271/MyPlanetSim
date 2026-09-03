#include "myplanetsim/numerics/compatible_operators.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "myplanetsim/numerics/spherical_operators.hpp"

namespace mps {
namespace {

void validate_cell_values(const CubedSphereGrid& grid,
                          const std::span<const Real> values) {
  if (values.size() != grid.cell_count()) {
    throw std::invalid_argument("cell field size does not match grid");
  }
  if (!std::ranges::all_of(values,
                           [](const Real value) { return std::isfinite(value); })) {
    throw std::invalid_argument("cell field contains a non-finite value");
  }
}

void validate_edge_values(const CubedSphereGrid& grid,
                          const std::span<const Real> values) {
  if (values.size() != grid.edge_count()) {
    throw std::invalid_argument("edge field size does not match grid");
  }
  if (!std::ranges::all_of(values,
                           [](const Real value) { return std::isfinite(value); })) {
    throw std::invalid_argument("edge field contains a non-finite value");
  }
}

void validate_vertex_values(const CubedSphereGrid& grid,
                            const std::span<const Real> values) {
  if (values.size() != grid.vertex_count()) {
    throw std::invalid_argument("vertex field size does not match grid");
  }
  if (!std::ranges::all_of(values,
                           [](const Real value) { return std::isfinite(value); })) {
    throw std::invalid_argument("vertex field contains a non-finite value");
  }
}

// Weights that interpolate a cell field to an edge along the dual edge.
struct EdgeInterpolationWeights {
  Real left;
  Real right;
};

[[nodiscard]] EdgeInterpolationWeights interpolation_weights(
    const DualEdgeGeometry& metric) {
  const Real span = metric.left_center_to_edge_m + metric.right_center_to_edge_m;
  return {.left = metric.right_center_to_edge_m / span,
          .right = metric.left_center_to_edge_m / span};
}

[[nodiscard]] Vec3 dual_edge_direction(const CubedSphereGrid& grid,
                                       const EdgeGeometry& edge) {
  return normalize(project_tangent(
      grid.cell(edge.right_cell).center - grid.cell(edge.left_cell).center,
      edge.center));
}

}  // namespace

std::vector<CompatibleEdgeMetric> compatible_edge_metrics(
    const CubedSphereGrid& grid, const CubedSphereDualTopology& dual) {
  std::vector<CompatibleEdgeMetric> metrics(grid.edge_count());
  for (const auto& edge : grid.edges()) {
    const auto basis = edge_tangent_basis(edge);
    const Vec3 direction = dual_edge_direction(grid, edge);
    metrics[edge.id] = {
        .hodge_ratio = edge.length_m / dual.edges()[edge.id].dual_length_m,
        .normal_projection = dot(direction, basis.normal),
        .tangent_projection = dot(direction, basis.tangent),
    };
  }
  return metrics;
}

std::vector<Real> edge_normal_velocity(const CubedSphereGrid& grid,
                                       const CubedSphereDualTopology& dual,
                                       const std::span<const Vec3> cell_velocity) {
  if (cell_velocity.size() != grid.cell_count()) {
    throw std::invalid_argument("vector field size does not match grid");
  }
  if (!std::ranges::all_of(cell_velocity, is_finite)) {
    throw std::invalid_argument("vector field contains a non-finite value");
  }
  std::vector<Real> normal_velocity(grid.edge_count());
  for (const auto& edge : grid.edges()) {
    const auto weights = interpolation_weights(dual.edges()[edge.id]);
    const Vec3 velocity =
        weights.left * cell_velocity[grid.cell_index(edge.left_cell)] +
        weights.right * cell_velocity[grid.cell_index(edge.right_cell)];
    normal_velocity[edge.id] = dot(velocity, edge_tangent_basis(edge).normal);
  }
  return normal_velocity;
}

std::vector<Vec3> perot_cell_velocity(
    const CubedSphereGrid& grid, const std::span<const Real> edge_normal_velocity) {
  validate_edge_values(grid, edge_normal_velocity);
  std::vector<Vec3> velocity(grid.cell_count());
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
    const Vec3 center = grid.cells()[cell].center;
    Vec3 accumulated{};
    for (const std::size_t edge_id : grid.cell_edges(grid.cell_id(cell))) {
      const auto& edge = grid.edge(edge_id);
      const Real outward =
          static_cast<Real>(grid.edge_sign_for_cell(edge_id, grid.cell_id(cell))) *
          edge_normal_velocity[edge_id];
      const Vec3 displacement =
          grid.radius_m() * project_tangent(edge.center - center, center);
      accumulated = accumulated + (edge.length_m * outward) * displacement;
    }
    velocity[cell] = project_tangent(accumulated / grid.cells()[cell].area_m2, center);
  }
  return velocity;
}

std::vector<Real> vertex_relative_vorticity(
    const CubedSphereGrid& grid, const CubedSphereDualTopology& dual,
    const std::span<const Real> edge_normal_velocity) {
  validate_edge_values(grid, edge_normal_velocity);
  const auto metrics = compatible_edge_metrics(grid, dual);
  const auto cell_velocity = perot_cell_velocity(grid, edge_normal_velocity);
  std::vector<Real> dual_velocity(grid.edge_count());
  for (const auto& edge : grid.edges()) {
    const auto weights = interpolation_weights(dual.edges()[edge.id]);
    const Vec3 reconstructed =
        weights.left * cell_velocity[grid.cell_index(edge.left_cell)] +
        weights.right * cell_velocity[grid.cell_index(edge.right_cell)];
    const Real tangential = dot(reconstructed, edge_tangent_basis(edge).tangent);
    dual_velocity[edge.id] =
        metrics[edge.id].normal_projection * edge_normal_velocity[edge.id] +
        metrics[edge.id].tangent_projection * tangential;
  }
  std::vector<Real> vorticity(grid.vertex_count());
  for (const auto& dual_vertex : dual.vertices()) {
    const std::size_t vertex = dual_vertex.primal_vertex;
    Real circulation = 0.0;
    for (const std::size_t edge : dual_vertex.incident_edges) {
      circulation += static_cast<Real>(dual.edge_vertex_incidence(edge, vertex)) *
                     dual_velocity[edge] * dual.edges()[edge].dual_length_m;
    }
    vorticity[vertex] = circulation / dual_vertex.area_m2;
  }
  return vorticity;
}

std::vector<Real> interpolate_cells_to_vertices(
    const CubedSphereGrid& grid, const CubedSphereDualTopology& dual,
    const std::span<const Real> cell_values) {
  validate_cell_values(grid, cell_values);
  std::vector<Real> vertex_values(grid.vertex_count());
  for (const auto& dual_vertex : dual.vertices()) {
    Real weighted = 0.0;
    for (const std::size_t cell : dual_vertex.incident_cells) {
      weighted += 0.25 * grid.cells()[cell].area_m2 * cell_values[cell];
    }
    vertex_values[dual_vertex.primal_vertex] = weighted / dual_vertex.area_m2;
  }
  return vertex_values;
}

std::vector<Real> vertex_potential_vorticity(
    const CubedSphereGrid& grid, const CubedSphereDualTopology& dual,
    const std::span<const Real> cell_depth,
    const std::span<const Real> edge_normal_velocity,
    const Vec3 rotation_vector_rad_s) {
  if (!is_finite(rotation_vector_rad_s)) {
    throw std::invalid_argument("rotation vector is non-finite");
  }
  const auto depth = interpolate_cells_to_vertices(grid, dual, cell_depth);
  auto potential_vorticity =
      vertex_relative_vorticity(grid, dual, edge_normal_velocity);
  for (std::size_t vertex = 0; vertex < grid.vertex_count(); ++vertex) {
    if (!(depth[vertex] > 0.0)) {
      throw std::invalid_argument("potential vorticity requires a positive depth");
    }
    potential_vorticity[vertex] =
        (potential_vorticity[vertex] +
         2.0 * dot(rotation_vector_rad_s, grid.vertices()[vertex].position)) /
        depth[vertex];
  }
  return potential_vorticity;
}

std::vector<Real> edge_potential_vorticity(
    const CubedSphereGrid& grid,
    const std::span<const Real> vertex_potential_vorticity) {
  validate_vertex_values(grid, vertex_potential_vorticity);
  std::vector<Real> edge_values(grid.edge_count());
  for (const auto& edge : grid.edges()) {
    edge_values[edge.id] = 0.5 * (vertex_potential_vorticity[edge.first_vertex] +
                                  vertex_potential_vorticity[edge.second_vertex]);
  }
  return edge_values;
}

std::vector<Real> cell_kinetic_energy(
    const CubedSphereGrid& grid, const std::span<const Real> edge_normal_velocity) {
  const auto velocity = perot_cell_velocity(grid, edge_normal_velocity);
  std::vector<Real> kinetic_energy(grid.cell_count());
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
    kinetic_energy[cell] = 0.5 * norm_squared(velocity[cell]);
  }
  return kinetic_energy;
}

}  // namespace mps
