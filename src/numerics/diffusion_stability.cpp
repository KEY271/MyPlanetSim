#include "myplanetsim/numerics/diffusion_stability.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "myplanetsim/core/validation.hpp"

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

}  // namespace mps
