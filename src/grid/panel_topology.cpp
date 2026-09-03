#include "myplanetsim/grid/panel_topology.hpp"

#include <cmath>
#include <stdexcept>

namespace mps {
namespace {

constexpr Real kQuarterPi = 0.785398163397448309615660845819875721;

[[nodiscard]] PanelCoordinates outside_point(const Panel panel, const PanelEdge edge,
                                             const Real parameter) {
  constexpr Real offset = 1.0e-7;
  Real alpha = parameter;
  Real beta = parameter;
  switch (edge) {
    case PanelEdge::kWest:
      alpha = -kQuarterPi - offset;
      break;
    case PanelEdge::kEast:
      alpha = kQuarterPi + offset;
      break;
    case PanelEdge::kSouth:
      beta = -kQuarterPi - offset;
      break;
    case PanelEdge::kNorth:
      beta = kQuarterPi + offset;
      break;
  }
  return inverse_map(map_to_unit_sphere(panel, alpha, beta));
}

[[nodiscard]] PanelEdge boundary_edge(const PanelCoordinates& coordinates) {
  const Real west = std::abs(coordinates.alpha + kQuarterPi);
  const Real east = std::abs(coordinates.alpha - kQuarterPi);
  const Real south = std::abs(coordinates.beta + kQuarterPi);
  const Real north = std::abs(coordinates.beta - kQuarterPi);
  if (west <= east && west <= south && west <= north) {
    return PanelEdge::kWest;
  }
  if (east <= south && east <= north) {
    return PanelEdge::kEast;
  }
  return south <= north ? PanelEdge::kSouth : PanelEdge::kNorth;
}

[[nodiscard]] Real edge_parameter(const PanelCoordinates& coordinates,
                                  const PanelEdge edge) {
  return edge == PanelEdge::kWest || edge == PanelEdge::kEast ? coordinates.beta
                                                              : coordinates.alpha;
}

}  // namespace

EdgeConnection edge_connection(const Panel panel, const PanelEdge edge) {
  const auto low = outside_point(panel, edge, -0.25);
  const auto high = outside_point(panel, edge, 0.25);
  if (low.panel != high.panel) {
    throw std::logic_error("panel edge does not map to one neighbor");
  }
  const auto neighbor_edge = boundary_edge(low);
  return {low.panel, neighbor_edge,
          edge_parameter(high, neighbor_edge) < edge_parameter(low, neighbor_edge)};
}

Index map_edge_index(const Index index, const Index cells_per_panel,
                     const EdgeConnection& connection) {
  if (cells_per_panel <= 0 || index < 0 || index >= cells_per_panel) {
    throw std::out_of_range("panel edge index is out of range");
  }
  return connection.reversed ? cells_per_panel - 1 - index : index;
}

}  // namespace mps
