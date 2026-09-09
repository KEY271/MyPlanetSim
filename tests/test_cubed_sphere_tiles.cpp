#include <algorithm>
#include <vector>

#include "myplanetsim/grid/cubed_sphere_tiles.hpp"
#include "support/test.hpp"

MPS_TEST_CASE("cubed-sphere tiles cover cells and own seam edges exactly once") {
  const mps::CubedSphereGrid grid(9, 2.0);
  for (const auto width : {8U, 16U, 32U}) {
    const auto tiles =
        mps::make_cubed_sphere_tiles(grid, width, mps::kDryHydrostaticTileHaloRadius);
    std::vector<unsigned int> cell_count(grid.cell_count(), 0);
    std::vector<unsigned int> edge_count(grid.edge_count(), 0);
    for (const auto& tile : tiles) {
      for (const auto cell : tile.interior_cells) ++cell_count[cell];
      for (const auto edge : tile.owned_edges) ++edge_count[edge];
      for (const auto cell : tile.interior_cells) {
        for (const auto& incident : grid.cell_cache()[cell].edges) {
          const auto& edge = grid.edge_cache()[incident.edge];
          const auto neighbor =
              edge.left_cell == cell ? edge.right_cell : edge.left_cell;
          const bool local = std::ranges::find(tile.interior_cells, neighbor) !=
                             tile.interior_cells.end();
          const bool halo =
              std::ranges::find(tile.halo_cells, neighbor) != tile.halo_cells.end();
          MPS_CHECK(local || halo);
        }
      }
    }
    for (const auto count : cell_count) MPS_CHECK_EQ(count, 1U);
    for (const auto count : edge_count) MPS_CHECK_EQ(count, 1U);
  }
}

MPS_TEST_CASE("tile halo crosses panel seams and reaches the requested radius") {
  const mps::CubedSphereGrid grid(9, 2.0);
  const auto tiles = mps::make_cubed_sphere_tiles(grid, 8, 2);
  const auto& corner_tile = tiles.front();
  MPS_CHECK(std::ranges::any_of(corner_tile.halo_cells, [&](const std::size_t cell) {
    return grid.cell_id(cell).panel != corner_tile.panel;
  }));
  MPS_CHECK(corner_tile.halo_cells.size() > 2 * corner_tile.interior_cells.size() / 9);
}

int main() { return mps::test::run_all(); }
