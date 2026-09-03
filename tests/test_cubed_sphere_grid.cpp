#include <cmath>
#include <numbers>

#include "myplanetsim/grid/cubed_sphere_grid.hpp"
#include "myplanetsim/grid/panel_topology.hpp"
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

MPS_TEST_CASE("panel edge connections are reciprocal") {
  for (const auto panel : mps::kPanels) {
    for (const auto edge : {mps::PanelEdge::kWest, mps::PanelEdge::kEast,
                            mps::PanelEdge::kSouth, mps::PanelEdge::kNorth}) {
      const auto first = mps::edge_connection(panel, edge);
      const auto second =
          mps::edge_connection(first.neighbor_panel, first.neighbor_edge);
      MPS_CHECK(second.neighbor_panel == panel);
      MPS_CHECK(second.neighbor_edge == edge);
      MPS_CHECK_EQ(second.reversed, first.reversed);
      for (mps::Index index = 0; index < 7; ++index) {
        const auto mapped = mps::map_edge_index(index, 7, first);
        MPS_CHECK_EQ(mps::map_edge_index(mapped, 7, second), index);
      }
    }
  }
}

int main() { return mps::test::run_all(); }
