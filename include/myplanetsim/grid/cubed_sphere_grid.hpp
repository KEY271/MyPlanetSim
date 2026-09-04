#pragma once

#include <array>
#include <cstddef>
#include <span>
#include <vector>

#include "myplanetsim/core/types.hpp"
#include "myplanetsim/geometry/cubed_sphere.hpp"
#include "myplanetsim/grid/cubed_sphere_cache.hpp"

namespace mps {

struct CellId {
  Panel panel;
  Index i;
  Index j;

  [[nodiscard]] friend constexpr bool operator==(const CellId&,
                                                 const CellId&) = default;
};

struct VertexGeometry {
  Vec3 position;
};

struct CellGeometry {
  CellId id;
  Vec3 center;
  Real area_m2;
  std::array<std::size_t, 4> vertices;
};

struct EdgeGeometry {
  std::size_t id;
  std::size_t first_vertex;
  std::size_t second_vertex;
  CellId left_cell;
  CellId right_cell;
  Vec3 center;
  Vec3 outward_normal_from_left;
  Real length_m;
};

class CubedSphereGrid {
 public:
  CubedSphereGrid(Index cells_per_panel, Real radius_m);

  [[nodiscard]] Index cells_per_panel() const noexcept { return n_; }
  [[nodiscard]] Real radius_m() const noexcept { return radius_m_; }
  [[nodiscard]] std::size_t cell_count() const noexcept { return cells_.size(); }
  [[nodiscard]] std::size_t edge_count() const noexcept { return edges_.size(); }
  [[nodiscard]] std::size_t vertex_count() const noexcept { return vertices_.size(); }
  [[nodiscard]] std::span<const CellGeometry> cells() const noexcept { return cells_; }
  [[nodiscard]] std::span<const EdgeGeometry> edges() const noexcept { return edges_; }
  [[nodiscard]] std::span<const VertexGeometry> vertices() const noexcept {
    return vertices_;
  }
  [[nodiscard]] std::span<const CachedCellGeometry> cell_cache() const noexcept {
    return cell_cache_;
  }
  [[nodiscard]] std::span<const CachedEdgeGeometry> edge_cache() const noexcept {
    return edge_cache_;
  }

  [[nodiscard]] std::size_t cell_index(CellId cell) const;
  [[nodiscard]] CellId cell_id(std::size_t flat_index) const;
  [[nodiscard]] const CellGeometry& cell(CellId id) const;
  [[nodiscard]] const EdgeGeometry& edge(std::size_t id) const;
  [[nodiscard]] std::array<std::size_t, 4> cell_edges(CellId id) const;
  [[nodiscard]] int edge_sign_for_cell(std::size_t edge_id, CellId cell) const;
  [[nodiscard]] CellId neighbor_across(std::size_t edge_id, CellId cell) const;
  [[nodiscard]] Real total_area_m2() const noexcept;

 private:
  Index n_;
  Real radius_m_;
  std::vector<VertexGeometry> vertices_;
  std::vector<CellGeometry> cells_;
  std::vector<EdgeGeometry> edges_;
  std::vector<std::array<std::size_t, 4>> cell_edges_;
  std::vector<CachedCellGeometry> cell_cache_;
  std::vector<CachedEdgeGeometry> edge_cache_;
};

}  // namespace mps
