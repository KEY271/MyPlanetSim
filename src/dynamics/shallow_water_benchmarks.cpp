#include "myplanetsim/dynamics/shallow_water_benchmarks.hpp"

#include <cmath>
#include <stdexcept>
#include <vector>

#include "myplanetsim/core/validation.hpp"

namespace mps {

Real linear_wave_frequency(const Real radius_m, const Real gravity_m_s2,
                           const Real mean_depth_m) {
  require_positive(radius_m, "linear-wave radius");
  require_positive(gravity_m_s2, "linear-wave gravity");
  require_positive(mean_depth_m, "linear-wave mean depth");
  return std::sqrt(2.0 * gravity_m_s2 * mean_depth_m) / radius_m;
}

ShallowWaterState make_linear_wave_state(const CubedSphereGrid& grid, const Real time_s,
                                         const Real gravity_m_s2,
                                         const Real mean_depth_m,
                                         const Real relative_amplitude) {
  require_non_negative(time_s, "linear-wave time");
  require_positive(relative_amplitude, "linear-wave relative amplitude");
  if (!(relative_amplitude < 1.0)) {
    throw std::invalid_argument("linear-wave amplitude must be less than one");
  }
  const Real frequency =
      linear_wave_frequency(grid.radius_m(), gravity_m_s2, mean_depth_m);
  const Real amplitude = relative_amplitude * mean_depth_m;
  const Real depth_phase = std::cos(frequency * time_s);
  const Real velocity_phase = std::sin(frequency * time_s);
  ShallowWaterState state{.time_s = time_s,
                          .step = 0,
                          .depth = std::vector<Real>(grid.cell_count()),
                          .momentum = std::vector<Vec3>(grid.cell_count())};
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
    const Vec3 position = grid.cells()[cell].center;
    state.depth[cell] = mean_depth_m + amplitude * position.x * depth_phase;
    const Vec3 velocity = -(gravity_m_s2 * amplitude / (frequency * grid.radius_m())) *
                          velocity_phase * project_tangent({1.0, 0.0, 0.0}, position);
    state.momentum[cell] = state.depth[cell] * velocity;
  }
  return state;
}

ShallowWaterState make_geostrophic_adjustment_state(const CubedSphereGrid& grid,
                                                    const Real time_s,
                                                    const Real mean_depth_m,
                                                    const Real relative_amplitude) {
  require_non_negative(time_s, "geostrophic-adjustment time");
  require_positive(mean_depth_m, "geostrophic-adjustment mean depth");
  require_positive(relative_amplitude, "geostrophic-adjustment relative amplitude");
  ShallowWaterState state{.time_s = time_s,
                          .step = 0,
                          .depth = std::vector<Real>(grid.cell_count()),
                          .momentum = std::vector<Vec3>(grid.cell_count())};
  constexpr Vec3 center{1.0, 0.0, 0.0};
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
    const Real angle = safe_angle(center, grid.cells()[cell].center);
    state.depth[cell] =
        mean_depth_m * (1.0 + relative_amplitude * std::exp(-20.0 * angle * angle));
  }
  return state;
}

ShallowWaterState make_williamson2_state(const CubedSphereGrid& grid, const Real time_s,
                                         const Real gravity_m_s2,
                                         const Real rotation_rate_rad_s,
                                         const Real reference_depth_m,
                                         const Real maximum_velocity_m_s,
                                         const Vec3 flow_axis) {
  require_non_negative(time_s, "Williamson 2 time");
  require_positive(gravity_m_s2, "Williamson 2 gravity");
  require_finite(rotation_rate_rad_s, "Williamson 2 rotation rate");
  require_positive(reference_depth_m, "Williamson 2 reference depth");
  require_positive(maximum_velocity_m_s, "Williamson 2 maximum velocity");
  if (!is_finite(flow_axis) || !(norm_squared(flow_axis) > 0.0)) {
    throw std::invalid_argument("Williamson 2 flow axis must be finite and nonzero");
  }
  const Vec3 axis = normalize(flow_axis);
  const Real height_coefficient =
      (grid.radius_m() * rotation_rate_rad_s * maximum_velocity_m_s +
       0.5 * maximum_velocity_m_s * maximum_velocity_m_s) /
      gravity_m_s2;
  if (!(reference_depth_m > height_coefficient)) {
    throw std::invalid_argument(
        "Williamson 2 reference depth does not keep depth positive");
  }
  ShallowWaterState state{.time_s = time_s,
                          .step = 0,
                          .depth = std::vector<Real>(grid.cell_count()),
                          .momentum = std::vector<Vec3>(grid.cell_count())};
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
    const Vec3 position = grid.cells()[cell].center;
    const Real sine_latitude = dot(axis, position);
    const Vec3 velocity = maximum_velocity_m_s * cross(axis, position);
    state.depth[cell] =
        reference_depth_m - height_coefficient * sine_latitude * sine_latitude;
    state.momentum[cell] = state.depth[cell] * velocity;
  }
  return state;
}

ShallowWaterState make_shallow_water_initial_state(const CubedSphereGrid& grid,
                                                   const ExperimentConfig& config) {
  switch (config.shallow_water.test_case) {
    case ShallowWaterTestCase::kRest:
      return {.time_s = config.run.start_time_s,
              .step = 0,
              .depth = std::vector<Real>(grid.cell_count(),
                                         config.shallow_water.mean_depth_m),
              .momentum = std::vector<Vec3>(grid.cell_count())};
    case ShallowWaterTestCase::kLinearWave:
      return make_linear_wave_state(grid, config.run.start_time_s,
                                    config.planet.gravity_m_s2,
                                    config.shallow_water.mean_depth_m);
    case ShallowWaterTestCase::kGeostrophicAdjustment:
      return make_geostrophic_adjustment_state(grid, config.run.start_time_s,
                                               config.shallow_water.mean_depth_m);
    case ShallowWaterTestCase::kWilliamson2:
      return make_williamson2_state(
          grid, config.run.start_time_s, config.planet.gravity_m_s2,
          config.planet.rotation_rate_rad_s, config.shallow_water.mean_depth_m,
          config.shallow_water.maximum_velocity_m_s,
          {config.shallow_water.flow_axis_x, config.shallow_water.flow_axis_y,
           config.shallow_water.flow_axis_z});
    case ShallowWaterTestCase::kWilliamson6:
    case ShallowWaterTestCase::kGalewsky:
      throw std::invalid_argument(
          "requested shallow-water benchmark is not yet available");
  }
  throw std::logic_error("unknown shallow-water test case");
}

}  // namespace mps
