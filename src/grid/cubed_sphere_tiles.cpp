#include "myplanetsim/grid/cubed_sphere_tiles.hpp"

#include <algorithm>
#include <limits>
#include <stdexcept>
#include <utility>

namespace mps {

std::vector<CubedSphereTile> make_cubed_sphere_tiles(const CubedSphereGrid& grid,
                                                     const std::size_t tile_width,
                                                     const std::size_t halo_radius) {
  if (tile_width == 0) throw std::invalid_argument("tile width must be positive");
  const auto n = static_cast<std::size_t>(grid.cells_per_panel());
  std::vector<CubedSphereTile> result;
  std::vector<std::size_t> cell_tile(grid.cell_count(),
                                     std::numeric_limits<std::size_t>::max());
  for (const Panel panel : kPanels) {
    for (std::size_t j = 0; j < n; j += tile_width) {
      for (std::size_t i = 0; i < n; i += tile_width) {
        auto& tile = result.emplace_back(
            CubedSphereTile{.id = result.size(),
                            .panel = panel,
                            .i_begin = static_cast<Index>(i),
                            .i_end = static_cast<Index>(std::min(n, i + tile_width)),
                            .j_begin = static_cast<Index>(j),
                            .j_end = static_cast<Index>(std::min(n, j + tile_width)),
                            .interior_cells = {},
                            .halo_cells = {},
                            .owned_edges = {}});
        for (Index cell_j = tile.j_begin; cell_j < tile.j_end; ++cell_j) {
          for (Index cell_i = tile.i_begin; cell_i < tile.i_end; ++cell_i) {
            const auto cell = grid.cell_index({panel, cell_i, cell_j});
            tile.interior_cells.push_back(cell);
            cell_tile[cell] = tile.id;
          }
        }
      }
    }
  }
  if (std::ranges::any_of(cell_tile, [](const std::size_t tile) {
        return tile == std::numeric_limits<std::size_t>::max();
      }))
    throw std::runtime_error("tile partition did not cover every cell");

  for (auto& tile : result) {
    std::vector<unsigned char> visited(grid.cell_count(), 0);
    std::vector<std::size_t> frontier = tile.interior_cells;
    for (const auto cell : frontier) visited[cell] = 1;
    for (std::size_t depth = 0; depth < halo_radius; ++depth) {
      std::vector<std::size_t> next;
      for (const auto cell : frontier) {
        for (const auto& incident : grid.cell_cache()[cell].edges) {
          const auto& edge = grid.edge_cache()[incident.edge];
          const auto neighbor =
              edge.left_cell == cell ? edge.right_cell : edge.left_cell;
          if (visited[neighbor] == 0) {
            visited[neighbor] = 1;
            next.push_back(neighbor);
          }
        }
      }
      frontier = std::move(next);
    }
    for (std::size_t cell = 0; cell < visited.size(); ++cell)
      if (visited[cell] != 0 && cell_tile[cell] != tile.id)
        tile.halo_cells.push_back(cell);
  }

  for (std::size_t edge_id = 0; edge_id < grid.edge_count(); ++edge_id) {
    const auto& edge = grid.edge_cache()[edge_id];
    const auto owner = std::min(cell_tile[edge.left_cell], cell_tile[edge.right_cell]);
    result[owner].owned_edges.push_back(edge_id);
  }
  return result;
}

}  // namespace mps
