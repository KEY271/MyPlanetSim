#pragma once

#include <span>
#include <vector>

#include "myplanetsim/geometry/vec3.hpp"
#include "myplanetsim/grid/cubed_sphere_grid.hpp"

namespace mps {

struct EdgeTangentBasis {
  Vec3 normal;
  Vec3 tangent;
};

struct TangentVectorGradient {
  Vec3 alpha_derivative;
  Vec3 beta_derivative;
};

[[nodiscard]] std::vector<Real> finite_volume_divergence(
    const CubedSphereGrid& grid, std::span<const Real> oriented_edge_flux);
[[nodiscard]] std::vector<Vec3> least_squares_gradient(
    const CubedSphereGrid& grid, std::span<const Real> cell_values);
[[nodiscard]] std::vector<Real> finite_volume_laplacian(
    const CubedSphereGrid& grid, std::span<const Real> cell_values);
[[nodiscard]] std::vector<Real> finite_volume_curl(const CubedSphereGrid& grid,
                                                   std::span<const Vec3> cell_vectors);
[[nodiscard]] std::vector<Real> finite_volume_vector_divergence(
    const CubedSphereGrid& grid, std::span<const Vec3> cell_vectors);
[[nodiscard]] std::vector<TangentVectorGradient> least_squares_vector_gradient(
    const CubedSphereGrid& grid, std::span<const Vec3> cell_vectors);
[[nodiscard]] Vec3 reconstruct_tangent_vector(const CubedSphereGrid& grid,
                                              std::size_t cell, Vec3 cell_value,
                                              Vec3 face_position,
                                              const TangentVectorGradient& gradient);
[[nodiscard]] Vec3 reconstruct_tangent_vector_cached(
    const CubedSphereGrid& grid, std::size_t cell, Vec3 cell_value, Vec3 face_position,
    Vec3 face_displacement_m, const TangentVectorGradient& gradient);
[[nodiscard]] EdgeTangentBasis edge_tangent_basis(const EdgeGeometry& edge);
[[nodiscard]] std::vector<Real> shallow_water_potential_vorticity(
    const CubedSphereGrid& grid, std::span<const Real> depth,
    std::span<const Vec3> velocity, Vec3 rotation_vector_rad_s);
[[nodiscard]] std::vector<Vec3> finite_volume_vector_laplacian(
    const CubedSphereGrid& grid, std::span<const Vec3> cell_vectors);
void scatter_oriented_edge_flux(const CubedSphereGrid& grid,
                                std::span<const Real> edge_flux,
                                std::span<Real> cell_tendency);

}  // namespace mps
