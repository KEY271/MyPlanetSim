#pragma once

#include "myplanetsim/core/types.hpp"
#include "myplanetsim/geometry/cubed_sphere.hpp"

namespace mps {

enum class PanelEdge { kWest, kEast, kSouth, kNorth };

struct EdgeConnection {
  Panel neighbor_panel;
  PanelEdge neighbor_edge;
  bool reversed;
};

[[nodiscard]] EdgeConnection edge_connection(Panel panel, PanelEdge edge);
[[nodiscard]] Index map_edge_index(Index index, Index cells_per_panel,
                                   const EdgeConnection& connection);

}  // namespace mps
