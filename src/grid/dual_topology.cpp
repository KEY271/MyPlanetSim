#include "myplanetsim/grid/dual_topology.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <utility>

namespace mps {
namespace {

[[nodiscard]] Vec3 ring_basis(const Vec3 position) {
  const Vec3 reference =
      std::abs(position.z) < 0.9 ? Vec3{0.0, 0.0, 1.0} : Vec3{1.0, 0.0, 0.0};
  return normalize(cross(reference, position));
}

template <typename Position>
void sort_around_vertex(std::vector<std::size_t>& ids, const Vec3 vertex,
                        Position&& position) {
  const Vec3 first_basis = ring_basis(vertex);
  const Vec3 second_basis = cross(vertex, first_basis);
  std::ranges::sort(ids, [&](const std::size_t first, const std::size_t second) {
    const Vec3 first_offset = project_tangent(position(first) - vertex, vertex);
    const Vec3 second_offset = project_tangent(position(second) - vertex, vertex);
    const Real first_angle =
        std::atan2(dot(first_offset, second_basis), dot(first_offset, first_basis));
    const Real second_angle =
        std::atan2(dot(second_offset, second_basis), dot(second_offset, first_basis));
    return first_angle < second_angle;
  });
}

}  // namespace

CubedSphereDualTopology::CubedSphereDualTopology(const CubedSphereGrid& grid)
    : grid_(&grid), vertices_(grid.vertex_count()), edges_(grid.edge_count()) {
  for (std::size_t vertex = 0; vertex < grid.vertex_count(); ++vertex) {
    vertices_[vertex].primal_vertex = vertex;
    vertices_[vertex].area_m2 = 0.0;
  }
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
    for (const std::size_t vertex : grid.cells()[cell].vertices) {
      vertices_[vertex].incident_cells.push_back(cell);
      vertices_[vertex].area_m2 += 0.25 * grid.cells()[cell].area_m2;
    }
  }
  for (const auto& edge : grid.edges()) {
    vertices_[edge.first_vertex].incident_edges.push_back(edge.id);
    vertices_[edge.second_vertex].incident_edges.push_back(edge.id);
    const Vec3 left = grid.cell(edge.left_cell).center;
    const Vec3 right = grid.cell(edge.right_cell).center;
    edges_[edge.id] = {
        .primal_edge = edge.id,
        .dual_length_m = grid.radius_m() * safe_angle(left, right),
        .left_center_to_edge_m = grid.radius_m() * safe_angle(left, edge.center),
        .right_center_to_edge_m = grid.radius_m() * safe_angle(right, edge.center),
    };
  }
  for (auto& dual_vertex : vertices_) {
    const Vec3 vertex = grid.vertices()[dual_vertex.primal_vertex].position;
    sort_around_vertex(dual_vertex.incident_cells, vertex, [&](const std::size_t cell) {
      return grid.cells()[cell].center;
    });
    sort_around_vertex(dual_vertex.incident_edges, vertex, [&](const std::size_t edge) {
      return grid.edges()[edge].center;
    });
    const std::size_t degree = dual_vertex.incident_cells.size();
    if ((degree != 3 && degree != 4) || dual_vertex.incident_edges.size() != degree ||
        !(dual_vertex.area_m2 > 0.0) || !std::isfinite(dual_vertex.area_m2)) {
      throw std::runtime_error("invalid cubed-sphere dual vertex");
    }
  }
  for (const auto& edge : edges_) {
    if (!(edge.dual_length_m > 0.0) || !(edge.left_center_to_edge_m > 0.0) ||
        !(edge.right_center_to_edge_m > 0.0) || !std::isfinite(edge.dual_length_m)) {
      throw std::runtime_error("invalid cubed-sphere dual edge metric");
    }
  }
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
    for (const std::size_t vertex : grid.cells()[cell].vertices) {
      if (boundary_of_boundary(cell, vertex) != 0) {
        throw std::runtime_error("cubed-sphere incidence violates D2 D1 = 0");
      }
    }
  }
}

int CubedSphereDualTopology::cell_edge_incidence(const std::size_t cell,
                                                 const std::size_t edge) const {
  if (cell >= grid_->cell_count()) {
    throw std::out_of_range("dual topology cell index is out of range");
  }
  const auto id = grid_->cell_id(cell);
  const auto& geometry = grid_->edge(edge);
  if (geometry.left_cell == id) {
    return 1;
  }
  if (geometry.right_cell == id) {
    return -1;
  }
  return 0;
}

int CubedSphereDualTopology::edge_vertex_incidence(const std::size_t edge,
                                                   const std::size_t vertex) const {
  if (vertex >= grid_->vertex_count()) {
    throw std::out_of_range("dual topology vertex index is out of range");
  }
  const auto& geometry = grid_->edge(edge);
  if (geometry.first_vertex == vertex) {
    return -1;
  }
  if (geometry.second_vertex == vertex) {
    return 1;
  }
  return 0;
}

int CubedSphereDualTopology::boundary_of_boundary(const std::size_t cell,
                                                  const std::size_t vertex) const {
  if (cell >= grid_->cell_count() || vertex >= grid_->vertex_count()) {
    throw std::out_of_range("dual topology incidence index is out of range");
  }
  int result = 0;
  for (const std::size_t edge : grid_->cell_edges(grid_->cell_id(cell))) {
    result += cell_edge_incidence(cell, edge) * edge_vertex_incidence(edge, vertex);
  }
  return result;
}

Real CubedSphereDualTopology::total_dual_area_m2() const noexcept {
  return std::accumulate(vertices_.begin(), vertices_.end(), Real{0.0},
                         [](const Real sum, const DualVertexGeometry& vertex) {
                           return sum + vertex.area_m2;
                         });
}

}  // namespace mps
