#include "myplanetsim/dynamics/shallow_water_state.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "myplanetsim/core/validation.hpp"

namespace mps {

Vec3 ShallowWaterState::velocity(const std::size_t cell) const {
  if (cell >= depth.size() || momentum.size() != depth.size()) {
    throw std::out_of_range("shallow-water cell index is out of range");
  }
  if (!(depth[cell] > 0.0) || !std::isfinite(depth[cell])) {
    throw std::runtime_error("cannot diagnose velocity from invalid depth");
  }
  return momentum[cell] / depth[cell];
}

std::vector<Real> flatten_shallow_water_state(const ShallowWaterState& state) {
  if (state.depth.empty() || state.momentum.size() != state.depth.size()) {
    throw std::invalid_argument("shallow-water state shape is invalid");
  }
  const std::size_t count = state.depth.size();
  std::vector<Real> flat(4 * count);
  std::copy(state.depth.begin(), state.depth.end(), flat.begin());
  for (std::size_t cell = 0; cell < count; ++cell) {
    flat[count + cell] = state.momentum[cell].x;
    flat[2 * count + cell] = state.momentum[cell].y;
    flat[3 * count + cell] = state.momentum[cell].z;
  }
  return flat;
}

ShallowWaterState unflatten_shallow_water_state(const Real time_s,
                                                const std::uint64_t step,
                                                const std::span<const Real> flat_state,
                                                const std::size_t cell_count) {
  if (cell_count == 0 || cell_count > std::numeric_limits<std::size_t>::max() / 4 ||
      flat_state.size() != 4 * cell_count) {
    throw std::invalid_argument("shallow-water flat state size does not match shape");
  }
  ShallowWaterState state{
      .time_s = time_s,
      .step = step,
      .depth = std::vector<Real>(flat_state.begin(), flat_state.begin() + cell_count),
      .momentum = std::vector<Vec3>(cell_count)};
  for (std::size_t cell = 0; cell < cell_count; ++cell) {
    state.momentum[cell] = {flat_state[cell_count + cell],
                            flat_state[2 * cell_count + cell],
                            flat_state[3 * cell_count + cell]};
  }
  return state;
}

void validate_shallow_water_state(const CubedSphereGrid& grid,
                                  const ShallowWaterState& state,
                                  const Real depth_floor_m, const Real momentum_scale) {
  require_non_negative(state.time_s, "shallow-water time_s");
  require_non_negative(depth_floor_m, "shallow-water depth floor");
  require_positive(momentum_scale, "shallow-water momentum scale");
  if (state.depth.size() != grid.cell_count() ||
      state.momentum.size() != grid.cell_count()) {
    throw std::invalid_argument("shallow-water state size does not match grid");
  }
  constexpr Real tolerance_factor = 64.0 * std::numeric_limits<Real>::epsilon();
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
    if (!std::isfinite(state.depth[cell]) || !(state.depth[cell] > depth_floor_m)) {
      throw std::runtime_error("shallow-water depth is non-finite or at its floor");
    }
    const Vec3 momentum = state.momentum[cell];
    if (!is_finite(momentum)) {
      throw std::runtime_error("shallow-water momentum is non-finite");
    }
    const Real scale = std::max(norm(momentum), momentum_scale);
    if (std::abs(dot(momentum, grid.cells()[cell].center)) > tolerance_factor * scale) {
      throw std::runtime_error("shallow-water momentum is not tangent");
    }
  }
}

void project_shallow_water_momentum(const CubedSphereGrid& grid,
                                    ShallowWaterState& state) {
  if (state.momentum.size() != grid.cell_count()) {
    throw std::invalid_argument("shallow-water momentum size does not match grid");
  }
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
    state.momentum[cell] =
        project_tangent(state.momentum[cell], grid.cells()[cell].center);
  }
}

}  // namespace mps
