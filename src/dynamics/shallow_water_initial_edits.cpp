#include "myplanetsim/dynamics/shallow_water_initial_edits.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <vector>

#include "myplanetsim/core/validation.hpp"

namespace mps {

void validate_initial_condition_edit(const InitialConditionEditV1& edit) {
  if (edit.kind != InitialConditionEditKind::kGaussianDepth) {
    throw std::invalid_argument("unsupported initial condition edit kind");
  }
  if (!is_finite(edit.center_unit)) {
    throw std::invalid_argument("initial edit center must be finite");
  }
  const Real center_norm = norm(edit.center_unit);
  if (!(center_norm > 0.0) || !std::isfinite(center_norm) ||
      std::abs(center_norm - 1.0) > 1.0e-12) {
    throw std::invalid_argument("initial edit center must be a unit vector");
  }
  require_finite(edit.amplitude_m, "initial edit amplitude_m");
  if (std::abs(edit.amplitude_m) > 1.0e6) {
    throw std::invalid_argument("initial edit amplitude_m is outside the supported range");
  }
  require_finite(edit.sigma_rad, "initial edit sigma_rad");
  if (!(edit.sigma_rad > 1.0e-6 && edit.sigma_rad <= std::numbers::pi_v<Real>)) {
    throw std::invalid_argument("initial edit sigma_rad is outside the supported range");
  }
}

InitialConditionEditDiagnostics apply_initial_condition_edits(
    const CubedSphereGrid& grid, ShallowWaterState& state,
    const std::span<const InitialConditionEditV1> edits, const Real depth_floor_m) {
  validate_shallow_water_state(grid, state, depth_floor_m);
  require_non_negative(depth_floor_m, "initial edit depth floor");
  Real total_added_volume = 0.0;
  for (const auto& edit : edits) {
    validate_initial_condition_edit(edit);
    std::vector<Real> perturbation(grid.cell_count());
    Real weighted_mean = 0.0;
    for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
      const Real cosine = std::clamp(dot(edit.center_unit, grid.cells()[cell].center), -1.0, 1.0);
      const Real distance = std::acos(cosine);
      perturbation[cell] = edit.amplitude_m *
                           std::exp(-(distance * distance) /
                                     (2.0 * edit.sigma_rad * edit.sigma_rad));
      weighted_mean += grid.cells()[cell].area_m2 * perturbation[cell];
    }
    if (edit.mass_policy == MassPolicy::kPreserveGlobal) {
      weighted_mean /= grid.total_area_m2();
    } else {
      weighted_mean = 0.0;
    }
    for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
      const Real change = perturbation[cell] - weighted_mean;
      state.depth[cell] += change;
      total_added_volume += grid.cells()[cell].area_m2 * change;
    }
    validate_shallow_water_state(grid, state, depth_floor_m);
  }
  const auto [minimum, maximum] = std::minmax_element(state.depth.begin(), state.depth.end());
  return {.added_volume_m3 = total_added_volume,
          .minimum_depth_m = *minimum,
          .maximum_depth_m = *maximum};
}

}  // namespace mps
