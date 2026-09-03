#include <algorithm>
#include <array>
#include <cmath>
#include <numbers>
#include <vector>

#include "myplanetsim/grid/cubed_sphere_grid.hpp"
#include "support/test.hpp"

MPS_TEST_CASE("cubed sphere geometry satisfies global invariants") {
  for (const mps::Index n : {1, 2, 4, 16, 64}) {
    constexpr double radius = 6.37122e6;
    const mps::CubedSphereGrid grid(n, radius);
    const auto count = static_cast<std::size_t>(n * n);
    MPS_CHECK_EQ(grid.cell_count(), 6U * count);
    MPS_CHECK_EQ(grid.edge_count(), 12U * count);
    MPS_CHECK_EQ(grid.vertex_count(), 6U * count + 2U);
    const double exact_area = 4.0 * std::numbers::pi * radius * radius;
    MPS_CHECK_NEAR(grid.total_area_m2(), exact_area, 5.0e-13 * exact_area);
    for (const auto& edge : grid.edges()) {
      MPS_CHECK(edge.length_m > 0.0);
      MPS_CHECK(edge.left_cell != edge.right_cell);
      MPS_CHECK_NEAR(mps::dot(edge.center, edge.outward_normal_from_left), 0.0,
                     5.0e-16);
      const auto& left = grid.cell(edge.left_cell);
      const auto& right = grid.cell(edge.right_cell);
      MPS_CHECK(mps::dot(edge.outward_normal_from_left, right.center - left.center) >
                0.0);
    }
  }
}

MPS_TEST_CASE("cubed sphere joins panel edges and all eight corners") {
  for (const mps::Index n : {1, 2, 7}) {
    const mps::CubedSphereGrid grid(n, 1.0);
    std::size_t panel_boundary_edges = 0;
    for (const auto& edge : grid.edges()) {
      if (edge.left_cell.panel == edge.right_cell.panel) {
        continue;
      }
      ++panel_boundary_edges;
      MPS_CHECK(grid.neighbor_across(edge.id, edge.left_cell) == edge.right_cell);
      MPS_CHECK(grid.neighbor_across(edge.id, edge.right_cell) == edge.left_cell);
      MPS_CHECK_EQ(grid.edge_sign_for_cell(edge.id, edge.left_cell), 1);
      MPS_CHECK_EQ(grid.edge_sign_for_cell(edge.id, edge.right_cell), -1);
    }
    MPS_CHECK_EQ(panel_boundary_edges, 12U * static_cast<std::size_t>(n));

    std::vector<std::array<bool, 6>> incident_panels(grid.vertex_count());
    for (const auto& cell : grid.cells()) {
      for (const auto vertex : cell.vertices) {
        incident_panels[vertex][mps::panel_index(cell.id.panel)] = true;
      }
    }
    std::size_t corners = 0;
    for (const auto& panels : incident_panels) {
      const auto count = std::ranges::count(panels, true);
      MPS_CHECK(count >= 1 && count <= 3);
      corners += count == 3 ? 1U : 0U;
    }
    MPS_CHECK_EQ(corners, 8U);
  }
}

int main() { return mps::test::run_all(); }
