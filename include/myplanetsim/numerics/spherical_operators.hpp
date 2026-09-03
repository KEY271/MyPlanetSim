#pragma once

#include <span>
#include <vector>

#include "myplanetsim/geometry/vec3.hpp"
#include "myplanetsim/grid/cubed_sphere_grid.hpp"

namespace mps {

[[nodiscard]] std::vector<Real> finite_volume_divergence(
    const CubedSphereGrid& grid, std::span<const Real> oriented_edge_flux);
[[nodiscard]] std::vector<Vec3> least_squares_gradient(
    const CubedSphereGrid& grid, std::span<const Real> cell_values);
[[nodiscard]] std::vector<Real> finite_volume_laplacian(
    const CubedSphereGrid& grid, std::span<const Real> cell_values);
[[nodiscard]] std::vector<Real> finite_volume_curl(const CubedSphereGrid& grid,
                                                   std::span<const Vec3> cell_vectors);
void scatter_oriented_edge_flux(const CubedSphereGrid& grid,
                                std::span<const Real> edge_flux,
                                std::span<Real> cell_tendency);

}  // namespace mps
