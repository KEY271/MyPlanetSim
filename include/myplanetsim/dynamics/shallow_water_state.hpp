#pragma once

#include <cstdint>
#include <span>
#include <string_view>
#include <vector>

#include "myplanetsim/geometry/vec3.hpp"
#include "myplanetsim/grid/cubed_sphere_grid.hpp"

namespace mps {

inline constexpr std::string_view kShallowWaterCheckpointLayout =
    "shallow_water_cell_v1";

struct ShallowWaterState {
  Real time_s = 0.0;
  std::uint64_t step = 0;
  std::vector<Real> depth;
  std::vector<Vec3> momentum;

  [[nodiscard]] std::size_t cell_count() const noexcept { return depth.size(); }
  [[nodiscard]] Vec3 velocity(std::size_t cell) const;
};

struct ShallowWaterTendency {
  std::vector<Real> depth;
  std::vector<Vec3> momentum;
};

[[nodiscard]] std::vector<Real> flatten_shallow_water_state(
    const ShallowWaterState& state);
[[nodiscard]] ShallowWaterState unflatten_shallow_water_state(
    Real time_s, std::uint64_t step, std::span<const Real> flat_state,
    std::size_t cell_count);
void validate_shallow_water_state(const CubedSphereGrid& grid,
                                  const ShallowWaterState& state, Real depth_floor_m,
                                  Real momentum_scale = 1.0);
void project_shallow_water_momentum(const CubedSphereGrid& grid,
                                    ShallowWaterState& state);

}  // namespace mps
