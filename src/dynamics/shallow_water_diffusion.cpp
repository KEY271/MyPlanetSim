#include "myplanetsim/dynamics/shallow_water_diffusion.hpp"

#include <stdexcept>
#include <vector>

#include "myplanetsim/core/validation.hpp"
#include "myplanetsim/numerics/spherical_operators.hpp"

namespace mps {
namespace {

[[nodiscard]] std::vector<Vec3> velocities(const CubedSphereGrid& grid,
                                           const ShallowWaterState& state) {
  std::vector<Vec3> result(grid.cell_count());
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
    result[cell] = state.velocity(cell);
  }
  return result;
}

}  // namespace

ShallowWaterTendency shallow_water_diffusion_tendency(
    const CubedSphereGrid& grid, const ShallowWaterState& state,
    const DiffusionKind kind, const Real diffusion_coefficient) {
  validate_shallow_water_state(grid, state, 0.0);
  require_non_negative(diffusion_coefficient, "diffusion coefficient");
  ShallowWaterTendency tendency{.depth = std::vector<Real>(grid.cell_count()),
                                .momentum = std::vector<Vec3>(grid.cell_count())};
  if (kind == DiffusionKind::kNone) {
    if (diffusion_coefficient != 0.0) {
      throw std::invalid_argument("disabled diffusion requires zero coefficient");
    }
    return tendency;
  }
  require_positive(diffusion_coefficient, "active diffusion coefficient");
  auto depth_operator = finite_volume_laplacian(grid, state.depth);
  auto velocity_operator =
      finite_volume_vector_laplacian(grid, velocities(grid, state));
  Real sign = 1.0;
  if (kind == DiffusionKind::kBiharmonic) {
    depth_operator = finite_volume_laplacian(grid, depth_operator);
    velocity_operator = finite_volume_vector_laplacian(grid, velocity_operator);
    sign = -1.0;
  }
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
    tendency.depth[cell] = sign * diffusion_coefficient * depth_operator[cell];
    const Vec3 velocity = state.velocity(cell);
    const Vec3 velocity_tendency =
        sign * diffusion_coefficient * velocity_operator[cell];
    tendency.momentum[cell] = project_tangent(
        state.depth[cell] * velocity_tendency + velocity * tendency.depth[cell],
        grid.cells()[cell].center);
  }
  return tendency;
}

}  // namespace mps
