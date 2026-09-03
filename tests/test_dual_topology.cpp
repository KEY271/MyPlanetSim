#include <algorithm>
#include <cmath>

#include "myplanetsim/grid/dual_topology.hpp"
#include "support/test.hpp"

MPS_TEST_CASE("cubed-sphere dual topology satisfies Euler and metric invariants") {
  for (const mps::Index n : {1, 2, 4, 16, 64}) {
    constexpr mps::Real radius = 6.37122e6;
    const mps::CubedSphereGrid grid(n, radius);
    const mps::CubedSphereDualTopology dual(grid);
    MPS_CHECK_EQ(grid.vertex_count() - grid.edge_count() + grid.cell_count(), 2U);
    MPS_CHECK_NEAR(dual.total_dual_area_m2(), grid.total_area_m2(),
                   5.0e-14 * grid.total_area_m2());
    std::size_t degree_three = 0;
    for (const auto& vertex : dual.vertices()) {
      MPS_CHECK(vertex.area_m2 > 0.0);
      MPS_CHECK_EQ(vertex.incident_cells.size(), vertex.incident_edges.size());
      MPS_CHECK(vertex.incident_cells.size() == 3U ||
                vertex.incident_cells.size() == 4U);
      degree_three += vertex.incident_cells.size() == 3U ? 1U : 0U;
    }
    MPS_CHECK_EQ(degree_three, 8U);
    for (const auto& edge : dual.edges()) {
      MPS_CHECK(edge.dual_length_m > 0.0);
      MPS_CHECK(edge.left_center_to_edge_m > 0.0);
      MPS_CHECK(edge.right_center_to_edge_m > 0.0);
    }
  }
}

MPS_TEST_CASE("signed cubed-sphere incidence has exact boundary of boundary") {
  for (const mps::Index n : {1, 2, 3, 8, 32}) {
    const mps::CubedSphereGrid grid(n, 1.0);
    const mps::CubedSphereDualTopology dual(grid);
    for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
      for (const std::size_t vertex : grid.cells()[cell].vertices) {
        MPS_CHECK_EQ(dual.boundary_of_boundary(cell, vertex), 0);
      }
      for (const std::size_t edge : grid.cell_edges(grid.cell_id(cell))) {
        const int incidence = dual.cell_edge_incidence(cell, edge);
        MPS_CHECK(incidence == -1 || incidence == 1);
      }
    }
    for (const auto& edge : grid.edges()) {
      MPS_CHECK_EQ(dual.edge_vertex_incidence(edge.id, edge.first_vertex), -1);
      MPS_CHECK_EQ(dual.edge_vertex_incidence(edge.id, edge.second_vertex), 1);
    }
  }
}

MPS_TEST_CASE("dual vertex rings are cyclically connected") {
  const mps::CubedSphereGrid grid(6, 1.0);
  const mps::CubedSphereDualTopology dual(grid);
  for (const auto& vertex : dual.vertices()) {
    for (std::size_t index = 0; index < vertex.incident_cells.size(); ++index) {
      const auto first = grid.cell_id(vertex.incident_cells[index]);
      const auto second = grid.cell_id(
          vertex.incident_cells[(index + 1) % vertex.incident_cells.size()]);
      bool share_edge = false;
      for (const std::size_t edge : grid.cell_edges(first)) {
        share_edge = share_edge || grid.neighbor_across(edge, first) == second;
      }
      MPS_CHECK(share_edge);
    }
  }
}

int main() { return mps::test::run_all(); }
