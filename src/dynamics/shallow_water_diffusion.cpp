#include "myplanetsim/dynamics/shallow_water_diffusion.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

#include "myplanetsim/core/validation.hpp"
#include "myplanetsim/numerics/spherical_operators.hpp"

namespace mps {
namespace {

[[nodiscard]] Real minimum_cell_spacing(const CubedSphereGrid& grid) {
  Real spacing = std::numeric_limits<Real>::infinity();
  for (const auto& edge : grid.edges()) {
    spacing = std::min(spacing,
                       grid.radius_m() * safe_angle(grid.cell(edge.left_cell).center,
                                                    grid.cell(edge.right_cell).center));
  }
  return spacing;
}

[[nodiscard]] std::vector<Vec3> velocities(const CubedSphereGrid& grid,
                                           const ShallowWaterState& state) {
  std::vector<Vec3> result(grid.cell_count());
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
    result[cell] = state.velocity(cell);
  }
  return result;
}

}  // namespace

Real stable_diffusion_time_step(const CubedSphereGrid& grid, const DiffusionKind kind,
                                const Real diffusion_coefficient,
                                const Real maximum_time_step_s) {
  require_positive(maximum_time_step_s, "maximum diffusion time step");
  require_non_negative(diffusion_coefficient, "diffusion coefficient");
  if (kind == DiffusionKind::kNone) {
    if (diffusion_coefficient != 0.0) {
      throw std::invalid_argument("disabled diffusion requires zero coefficient");
    }
    return maximum_time_step_s;
  }
  require_positive(diffusion_coefficient, "active diffusion coefficient");
  const Real spacing = minimum_cell_spacing(grid);
  const Real stability_limit =
      kind == DiffusionKind::kLaplacian
          ? 0.25 * spacing * spacing / diffusion_coefficient
          : 0.125 * std::pow(spacing, 4) / diffusion_coefficient;
  return std::min(maximum_time_step_s, stability_limit);
}

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
