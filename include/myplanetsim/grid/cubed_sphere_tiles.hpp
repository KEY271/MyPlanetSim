#pragma once

#include <cstddef>
#include <vector>

#include "myplanetsim/grid/cubed_sphere_grid.hpp"

namespace mps {

struct CubedSphereTile {
  std::size_t id = 0;
  Panel panel = Panel::kPositiveX;
  Index i_begin = 0;
  Index i_end = 0;
  Index j_begin = 0;
  Index j_end = 0;
  std::vector<std::size_t> interior_cells;
  std::vector<std::size_t> halo_cells;
  std::vector<std::size_t> owned_edges;
};

// Halo radius two covers the widest currently implemented horizontal stencil:
// the two-stage Laplacian used by explicit diffusion.
inline constexpr std::size_t kDryHydrostaticTileHaloRadius = 2;
inline constexpr std::size_t kDryHydrostaticTileWidth = 16;

[[nodiscard]] std::vector<CubedSphereTile> make_cubed_sphere_tiles(
    const CubedSphereGrid& grid, std::size_t tile_width, std::size_t halo_radius);

}  // namespace mps
