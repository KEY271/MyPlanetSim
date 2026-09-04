#pragma once

#include <array>
#include <cstddef>

#include "myplanetsim/geometry/cubed_sphere.hpp"

namespace mps {

struct CachedCellEdgeGeometry {
  std::size_t edge;
  std::size_t neighbor;
  int sign;
  Vec3 outward_normal;
  Vec3 neighbor_displacement_m;
  Vec3 face_displacement_m;
  Vec3 normalized_face_displacement_m;
  TangentComponents neighbor_coordinates_m;
  TangentComponents face_coordinates_m;
};

struct CachedCellGeometry {
  TangentBasis basis;
  Real inverse_area_m2;
  std::array<CachedCellEdgeGeometry, 4> edges;
  Vec3 pressure_geometry_correction_m;
};

struct CachedEdgeGeometry {
  std::size_t left_cell;
  std::size_t right_cell;
  Vec3 normal;
  Vec3 tangent;
  Vec3 circulation_tangent;
  Real center_distance_m;
};

}  // namespace mps
