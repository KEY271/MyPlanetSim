#include "myplanetsim/grid/cubed_sphere_grid.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <numbers>
#include <stdexcept>
#include <tuple>
#include <utility>

#include "myplanetsim/core/validation.hpp"

namespace mps {
namespace {

using VertexKey = std::tuple<std::int64_t, std::int64_t, std::int64_t>;
using EdgeKey = std::pair<std::size_t, std::size_t>;

[[nodiscard]] VertexKey vertex_key(const Vec3 position) {
  constexpr Real scale = 1.0e12;
  return {static_cast<std::int64_t>(std::floor(position.x * scale)),
          static_cast<std::int64_t>(std::floor(position.y * scale)),
          static_cast<std::int64_t>(std::floor(position.z * scale))};
}

[[nodiscard]] Real triangle_solid_angle(const Vec3 a, const Vec3 b,
                                        const Vec3 c) noexcept {
  return 2.0 * std::atan2(dot(a, cross(b, c)), 1.0 + dot(a, b) + dot(b, c) + dot(c, a));
}

[[nodiscard]] Real cell_area(const std::array<Vec3, 4>& corners,
                             const Real radius_m) noexcept {
  const Real solid_angle = triangle_solid_angle(corners[0], corners[1], corners[2]) +
                           triangle_solid_angle(corners[0], corners[2], corners[3]);
  return std::abs(solid_angle) * radius_m * radius_m;
}

[[nodiscard]] Vec3 logarithmic_displacement(const Vec3 from, const Vec3 to,
                                            const Real radius_m) {
  const Real cosine = std::clamp(dot(from, to), -1.0, 1.0);
  const Real angle = std::acos(cosine);
  const Real sine = std::sin(angle);
  if (!(sine > 0.0)) {
    throw std::runtime_error("degenerate cached grid displacement");
  }
  return (radius_m * angle / sine) * (to - cosine * from);
}

}  // namespace

CubedSphereGrid::CubedSphereGrid(const Index cells_per_panel, const Real radius_m)
    : n_(cells_per_panel), radius_m_(radius_m) {
  if (n_ <= 0) {
    throw std::invalid_argument("cells_per_panel must be positive");
  }
  require_positive(radius_m_, "radius_m");
  const auto n = static_cast<std::size_t>(n_);
  if (n > std::numeric_limits<std::size_t>::max() / n ||
      n * n > std::numeric_limits<std::size_t>::max() / 12) {
    throw std::length_error("cubed-sphere grid size overflows size_t");
  }

  const Real delta = std::numbers::pi_v<Real> / (2.0 * static_cast<Real>(n_));
  std::map<VertexKey, std::size_t> vertex_ids;
  std::array<std::vector<std::size_t>, 6> panel_vertices;
  for (const Panel panel : kPanels) {
    auto& ids = panel_vertices[panel_index(panel)];
    ids.resize((n + 1) * (n + 1));
    for (Index j = 0; j <= n_; ++j) {
      for (Index i = 0; i <= n_; ++i) {
        const Real alpha =
            -std::numbers::pi_v<Real> / 4.0 + static_cast<Real>(i) * delta;
        const Real beta =
            -std::numbers::pi_v<Real> / 4.0 + static_cast<Real>(j) * delta;
        const Vec3 position = map_to_unit_sphere(panel, alpha, beta);
        const auto [key_x, key_y, key_z] = vertex_key(position);
        std::size_t vertex_id = vertices_.size();
        for (std::int64_t dz = -1; dz <= 1 && vertex_id == vertices_.size(); ++dz) {
          for (std::int64_t dy = -1; dy <= 1 && vertex_id == vertices_.size(); ++dy) {
            for (std::int64_t dx = -1; dx <= 1; ++dx) {
              const auto found = vertex_ids.find({key_x + dx, key_y + dy, key_z + dz});
              if (found != vertex_ids.end() &&
                  norm(vertices_[found->second].position - position) < 1.0e-11) {
                vertex_id = found->second;
                break;
              }
            }
          }
        }
        if (vertex_id == vertices_.size()) {
          vertex_id = vertices_.size();
          const auto [iterator, inserted] =
              vertex_ids.emplace(VertexKey{key_x, key_y, key_z}, vertex_id);
          if (!inserted) {
            throw std::runtime_error("cubed-sphere vertex key collision");
          }
          static_cast<void>(iterator);
          vertices_.push_back({position});
        }
        ids[static_cast<std::size_t>(j) * (n + 1) + static_cast<std::size_t>(i)] =
            vertex_id;
      }
    }
  }

  cells_.reserve(6 * n * n);
  cell_edges_.resize(6 * n * n);
  std::map<EdgeKey, std::size_t> edge_ids;
  for (const Panel panel : kPanels) {
    const auto& ids = panel_vertices[panel_index(panel)];
    for (Index j = 0; j < n_; ++j) {
      for (Index i = 0; i < n_; ++i) {
        const auto vertex_at = [&](const Index vertex_i, const Index vertex_j) {
          return ids[static_cast<std::size_t>(vertex_j) * (n + 1) +
                     static_cast<std::size_t>(vertex_i)];
        };
        const std::array<std::size_t, 4> corner_ids{
            vertex_at(i, j), vertex_at(i + 1, j), vertex_at(i + 1, j + 1),
            vertex_at(i, j + 1)};
        std::array<Vec3, 4> corners{};
        for (std::size_t corner = 0; corner < 4; ++corner) {
          corners[corner] = vertices_[corner_ids[corner]].position;
        }
        const CellId id{panel, i, j};
        const Real alpha =
            -std::numbers::pi_v<Real> / 4.0 + (static_cast<Real>(i) + 0.5) * delta;
        const Real beta =
            -std::numbers::pi_v<Real> / 4.0 + (static_cast<Real>(j) + 0.5) * delta;
        cells_.push_back({id, map_to_unit_sphere(panel, alpha, beta),
                          cell_area(corners, radius_m_), corner_ids});

        const std::array<std::pair<std::size_t, std::size_t>, 4> directed_edges{
            {{corner_ids[0], corner_ids[1]},
             {corner_ids[1], corner_ids[2]},
             {corner_ids[2], corner_ids[3]},
             {corner_ids[3], corner_ids[0]}}};
        auto& local_edges = cell_edges_[cell_index(id)];
        for (std::size_t side = 0; side < 4; ++side) {
          const auto [first, second] = directed_edges[side];
          const EdgeKey key = std::minmax(first, second);
          const auto found = edge_ids.find(key);
          if (found == edge_ids.end()) {
            const std::size_t edge_id = edges_.size();
            edge_ids.emplace(key, edge_id);
            const Vec3 first_position = vertices_[first].position;
            const Vec3 second_position = vertices_[second].position;
            const Vec3 center = normalize(first_position + second_position);
            const Vec3 tangent =
                normalize(project_tangent(second_position - first_position, center));
            const Vec3 outward = normalize(cross(tangent, center));
            edges_.push_back({edge_id, first, second, id, id, center, outward,
                              radius_m_ * safe_angle(first_position, second_position)});
            local_edges[side] = edge_id;
          } else {
            auto& edge = edges_[found->second];
            if (edge.right_cell != edge.left_cell) {
              throw std::runtime_error("cubed-sphere edge has more than two cells");
            }
            if (edge.first_vertex != second || edge.second_vertex != first) {
              throw std::runtime_error("neighbor cells disagree on edge orientation");
            }
            edge.right_cell = id;
            local_edges[side] = found->second;
          }
        }
      }
    }
  }

  if (vertices_.size() != 6 * n * n + 2 || edges_.size() != 12 * n * n) {
    throw std::runtime_error("cubed-sphere topology count invariant failed");
  }
  for (const auto& cell_geometry : cells_) {
    if (!(cell_geometry.area_m2 > 0.0) || !std::isfinite(cell_geometry.area_m2)) {
      throw std::runtime_error("cubed-sphere cell has invalid area");
    }
  }
  for (const auto& edge_geometry : edges_) {
    if (edge_geometry.left_cell == edge_geometry.right_cell ||
        !(edge_geometry.length_m > 0.0) || !std::isfinite(edge_geometry.length_m)) {
      throw std::runtime_error("cubed-sphere edge is incomplete or invalid");
    }
  }

  edge_cache_.resize(edges_.size());
  for (const auto& edge_geometry : edges_) {
    const std::size_t left = cell_index(edge_geometry.left_cell);
    const std::size_t right = cell_index(edge_geometry.right_cell);
    const Vec3 normal = normalize(
        project_tangent(edge_geometry.outward_normal_from_left, edge_geometry.center));
    edge_cache_[edge_geometry.id] = {
        .left_cell = left,
        .right_cell = right,
        .normal = normal,
        .tangent = normalize(cross(edge_geometry.center, normal)),
        .circulation_tangent = normalize(
            project_tangent(vertices_[edge_geometry.second_vertex].position -
                                vertices_[edge_geometry.first_vertex].position,
                            edge_geometry.center)),
        .center_distance_m =
            radius_m_ * safe_angle(cells_[left].center, cells_[right].center),
    };
  }

  cell_cache_.resize(cells_.size());
  for (std::size_t cell = 0; cell < cells_.size(); ++cell) {
    const auto& geometry = cells_[cell];
    const auto coordinates = inverse_map(geometry.center);
    auto& cache = cell_cache_[cell];
    cache.basis = tangent_basis(coordinates.panel, coordinates.alpha, coordinates.beta);
    cache.inverse_area_m2 = 1.0 / geometry.area_m2;
    cache.pressure_geometry_correction_m = {};
    for (std::size_t side = 0; side < cache.edges.size(); ++side) {
      const std::size_t edge_id = cell_edges_[cell][side];
      const auto& edge_geometry = edges_[edge_id];
      const auto& cached_edge = edge_cache_[edge_id];
      const int sign = edge_geometry.left_cell == geometry.id ? 1 : -1;
      const std::size_t neighbor =
          sign == 1 ? cached_edge.right_cell : cached_edge.left_cell;
      const Vec3 neighbor_displacement =
          logarithmic_displacement(geometry.center, cells_[neighbor].center, radius_m_);
      const Vec3 face_displacement =
          logarithmic_displacement(geometry.center, edge_geometry.center, radius_m_);
      const Vec3 normalized_face_displacement = logarithmic_displacement(
          geometry.center, normalize(edge_geometry.center), radius_m_);
      cache.edges[side] = {
          .edge = edge_id,
          .neighbor = neighbor,
          .sign = sign,
          .outward_normal = normalize(project_tangent(
              static_cast<Real>(sign) * edge_geometry.outward_normal_from_left,
              geometry.center)),
          .neighbor_displacement_m = neighbor_displacement,
          .face_displacement_m = face_displacement,
          .normalized_face_displacement_m = normalized_face_displacement,
          .neighbor_coordinates_m = {dot(neighbor_displacement, cache.basis.alpha),
                                     dot(neighbor_displacement, cache.basis.beta)},
          .face_coordinates_m = {dot(face_displacement, cache.basis.alpha),
                                 dot(face_displacement, cache.basis.beta)},
      };
      if (sign == 1) {
        edge_cache_[edge_id].left_slot = side;
      } else {
        edge_cache_[edge_id].right_slot = side;
      }
      cache.pressure_geometry_correction_m =
          cache.pressure_geometry_correction_m +
          static_cast<Real>(sign) * edge_geometry.length_m * cached_edge.normal;
    }

    Real a00 = 0.0;
    Real a01 = 0.0;
    Real a11 = 0.0;
    for (const auto& edge : cache.edges) {
      const Real x = edge.neighbor_coordinates_m.alpha;
      const Real y = edge.neighbor_coordinates_m.beta;
      const Real weight = 1.0 / std::max(x * x + y * y, 1.0e-300);
      a00 += weight * x * x;
      a01 += weight * x * y;
      a11 += weight * y * y;
    }
    const Real determinant = a00 * a11 - a01 * a01;
    if (!(determinant > 1.0e-12 * a00 * a11)) {
      throw std::runtime_error("least-squares grid stencil is rank deficient");
    }
    for (auto& edge : cache.edges) {
      const Real x = edge.neighbor_coordinates_m.alpha;
      const Real y = edge.neighbor_coordinates_m.beta;
      const Real weight = 1.0 / std::max(x * x + y * y, 1.0e-300);
      edge.least_squares_weight_m_inverse = {
          .alpha = (a11 * weight * x - a01 * weight * y) / determinant,
          .beta = (a00 * weight * y - a01 * weight * x) / determinant,
      };
    }
  }
}

std::size_t CubedSphereGrid::cell_index(const CellId cell) const {
  if (cell.i < 0 || cell.i >= n_ || cell.j < 0 || cell.j >= n_ ||
      panel_index(cell.panel) >= kPanels.size()) {
    throw std::out_of_range("cubed-sphere cell index is out of range");
  }
  const auto n = static_cast<std::size_t>(n_);
  return (panel_index(cell.panel) * n + static_cast<std::size_t>(cell.j)) * n +
         static_cast<std::size_t>(cell.i);
}

CellId CubedSphereGrid::cell_id(const std::size_t flat_index) const {
  if (flat_index >= cells_.size()) {
    throw std::out_of_range("cubed-sphere flat cell index is out of range");
  }
  return cells_[flat_index].id;
}

const CellGeometry& CubedSphereGrid::cell(const CellId id) const {
  return cells_[cell_index(id)];
}

const EdgeGeometry& CubedSphereGrid::edge(const std::size_t id) const {
  if (id >= edges_.size()) {
    throw std::out_of_range("cubed-sphere edge index is out of range");
  }
  return edges_[id];
}

std::array<std::size_t, 4> CubedSphereGrid::cell_edges(const CellId id) const {
  return cell_edges_[cell_index(id)];
}

int CubedSphereGrid::edge_sign_for_cell(const std::size_t edge_id,
                                        const CellId cell_id_value) const {
  const auto& edge_geometry = edge(edge_id);
  if (edge_geometry.left_cell == cell_id_value) {
    return 1;
  }
  if (edge_geometry.right_cell == cell_id_value) {
    return -1;
  }
  throw std::invalid_argument("cell is not incident to edge");
}

CellId CubedSphereGrid::neighbor_across(const std::size_t edge_id,
                                        const CellId cell_id_value) const {
  const auto& edge_geometry = edge(edge_id);
  if (edge_geometry.left_cell == cell_id_value) {
    return edge_geometry.right_cell;
  }
  if (edge_geometry.right_cell == cell_id_value) {
    return edge_geometry.left_cell;
  }
  throw std::invalid_argument("cell is not incident to edge");
}

Real CubedSphereGrid::total_area_m2() const noexcept {
  Real sum = 0.0;
  for (const auto& cell_geometry : cells_) {
    sum += cell_geometry.area_m2;
  }
  return sum;
}

}  // namespace mps
