#pragma once

#include <cstddef>
#include <span>
#include <vector>

#include "myplanetsim/grid/cubed_sphere_grid.hpp"

namespace mps {

struct DualVertexGeometry {
  std::size_t primal_vertex;
  std::vector<std::size_t> incident_cells;
  std::vector<std::size_t> incident_edges;
  Real area_m2;
};

struct DualEdgeGeometry {
  std::size_t primal_edge;
  Real dual_length_m;
  Real left_center_to_edge_m;
  Real right_center_to_edge_m;
};

class CubedSphereDualTopology {
 public:
  explicit CubedSphereDualTopology(const CubedSphereGrid& grid);

  [[nodiscard]] std::span<const DualVertexGeometry> vertices() const noexcept {
    return vertices_;
  }
  [[nodiscard]] std::span<const DualEdgeGeometry> edges() const noexcept {
    return edges_;
  }
  [[nodiscard]] int cell_edge_incidence(std::size_t cell, std::size_t edge) const;
  [[nodiscard]] int edge_vertex_incidence(std::size_t edge, std::size_t vertex) const;
  [[nodiscard]] int boundary_of_boundary(std::size_t cell, std::size_t vertex) const;
  [[nodiscard]] Real total_dual_area_m2() const noexcept;

 private:
  const CubedSphereGrid* grid_;
  std::vector<DualVertexGeometry> vertices_;
  std::vector<DualEdgeGeometry> edges_;
};

}  // namespace mps
