#include "myplanetsim/dynamics/shallow_water_benchmarks.hpp"

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstdint>
#include <numbers>
#include <stdexcept>
#include <vector>

#include "myplanetsim/core/validation.hpp"

namespace mps {
namespace {

struct AxisBasis {
  Vec3 axis;
  Vec3 longitude_zero;
  Vec3 longitude_quarter;
};

[[nodiscard]] AxisBasis make_axis_basis(const Vec3 requested_axis) {
  if (!is_finite(requested_axis) || !(norm_squared(requested_axis) > 0.0)) {
    throw std::invalid_argument("benchmark axis must be finite and nonzero");
  }
  const Vec3 axis = normalize(requested_axis);
  const Vec3 reference =
      std::abs(axis.z) > 0.9 ? Vec3{1.0, 0.0, 0.0} : Vec3{0.0, 0.0, 1.0};
  const Vec3 longitude_zero = normalize(project_tangent(reference, axis));
  return {.axis = axis,
          .longitude_zero = longitude_zero,
          .longitude_quarter = cross(axis, longitude_zero)};
}

[[nodiscard]] Real longitude(const AxisBasis& basis, const Vec3 position) {
  return std::atan2(dot(position, basis.longitude_quarter),
                    dot(position, basis.longitude_zero));
}

void checksum_value(std::uint64_t& hash, const Real value) noexcept {
  constexpr std::uint64_t kPrime = 1099511628211ULL;
  const std::uint64_t bits = std::bit_cast<std::uint64_t>(value);
  for (int byte = 0; byte < 8; ++byte) {
    hash ^= (bits >> (8 * byte)) & 0xffU;
    hash *= kPrime;
  }
}

}  // namespace

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

ShallowWaterState make_williamson6_state(const CubedSphereGrid& grid, const Real time_s,
                                         const Real gravity_m_s2,
                                         const Real rotation_rate_rad_s,
                                         const Real reference_depth_m,
                                         const Real velocity_scale_m_s,
                                         const Vec3 flow_axis) {
  require_non_negative(time_s, "Williamson 6 time");
  require_positive(gravity_m_s2, "Williamson 6 gravity");
  require_finite(rotation_rate_rad_s, "Williamson 6 rotation rate");
  require_positive(reference_depth_m, "Williamson 6 reference depth");
  require_positive(velocity_scale_m_s, "Williamson 6 velocity scale");
  const AxisBasis basis = make_axis_basis(flow_axis);
  constexpr int wave_number = 4;
  const Real radius = grid.radius_m();
  const Real omega = velocity_scale_m_s / radius;
  const Real wave_rate = omega;
  ShallowWaterState state{.time_s = time_s,
                          .step = 0,
                          .depth = std::vector<Real>(grid.cell_count()),
                          .momentum = std::vector<Vec3>(grid.cell_count())};
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
    const Vec3 position = grid.cells()[cell].center;
    const Real sine_latitude = dot(basis.axis, position);
    const Real cosine_latitude =
        std::sqrt(std::max(0.0, 1.0 - sine_latitude * sine_latitude));
    const Real lambda = longitude(basis, position);
    const Real cosine_wave = std::cos(wave_number * lambda);
    const Real sine_wave = std::sin(wave_number * lambda);
    const Real cosine_squared = cosine_latitude * cosine_latitude;
    const Real cosine_power = std::pow(cosine_latitude, wave_number - 1);
    const Real zonal =
        radius * omega * cosine_latitude +
        radius * wave_rate * cosine_power *
            (wave_number * sine_latitude * sine_latitude - cosine_squared) *
            cosine_wave;
    const Real meridional =
        -radius * wave_rate * wave_number * cosine_power * sine_latitude * sine_wave;
    const Vec3 east = cosine_latitude > 1.0e-14
                          ? cross(basis.axis, position) / cosine_latitude
                          : basis.longitude_quarter;
    const Vec3 north =
        project_tangent(basis.axis, position) / std::max(cosine_latitude, 1.0e-14);

    const Real cosine_2r = std::pow(cosine_latitude, 2 * wave_number);
    const Real cosine_2r_minus_2 = std::pow(cosine_latitude, 2 * wave_number - 2);
    const Real a =
        0.5 * omega * (2.0 * rotation_rate_rad_s + omega) * cosine_squared +
        0.25 * wave_rate * wave_rate *
            (cosine_2r * ((wave_number + 1.0) * cosine_squared +
                          (2.0 * wave_number * wave_number - wave_number - 2.0)) -
             2.0 * wave_number * wave_number * cosine_2r_minus_2);
    const Real b = 2.0 * (rotation_rate_rad_s + omega) * wave_rate /
                   ((wave_number + 1.0) * (wave_number + 2.0)) *
                   std::pow(cosine_latitude, wave_number) *
                   ((wave_number * wave_number + 2.0 * wave_number + 2.0) -
                    (wave_number + 1.0) * (wave_number + 1.0) * cosine_squared);
    const Real c = 0.25 * wave_rate * wave_rate * cosine_2r *
                   ((wave_number + 1.0) * cosine_squared - (wave_number + 2.0));
    state.depth[cell] =
        reference_depth_m +
        radius * radius / gravity_m_s2 *
            (a + b * cosine_wave + c * std::cos(2.0 * wave_number * lambda));
    if (!(state.depth[cell] > 0.0) || !std::isfinite(state.depth[cell])) {
      throw std::invalid_argument("Williamson 6 constants produce invalid depth");
    }
    state.momentum[cell] = state.depth[cell] * (zonal * east + meridional * north);
  }
  return state;
}

ShallowWaterState make_galewsky_state(const CubedSphereGrid& grid, const Real time_s,
                                      const Real gravity_m_s2,
                                      const Real rotation_rate_rad_s,
                                      const Real mean_depth_m,
                                      const Real maximum_velocity_m_s,
                                      const Vec3 flow_axis,
                                      const bool add_height_perturbation) {
  require_non_negative(time_s, "Galewsky time");
  require_positive(gravity_m_s2, "Galewsky gravity");
  require_finite(rotation_rate_rad_s, "Galewsky rotation rate");
  require_positive(mean_depth_m, "Galewsky mean depth");
  require_positive(maximum_velocity_m_s, "Galewsky maximum velocity");
  const AxisBasis basis = make_axis_basis(flow_axis);
  constexpr Real phi0 = std::numbers::pi_v<Real> / 7.0;
  constexpr Real phi1 = 0.5 * std::numbers::pi_v<Real> - phi0;
  constexpr Real perturbation_latitude = std::numbers::pi_v<Real> / 4.0;
  constexpr Real perturbation_longitude_scale = 1.0 / 3.0;
  constexpr Real perturbation_latitude_scale = 1.0 / 15.0;
  constexpr Real perturbation_height_m = 120.0;
  constexpr std::size_t integration_intervals = 4096;
  constexpr Real latitude_min = -0.5 * std::numbers::pi_v<Real>;
  constexpr Real latitude_step =
      std::numbers::pi_v<Real> / static_cast<Real>(integration_intervals);
  const Real normalization = std::exp(-4.0 / ((phi1 - phi0) * (phi1 - phi0)));
  const auto zonal_velocity = [&](const Real latitude) {
    if (!(latitude > phi0 && latitude < phi1)) {
      return 0.0;
    }
    return maximum_velocity_m_s / normalization *
           std::exp(1.0 / ((latitude - phi0) * (latitude - phi1)));
  };
  const auto balance_integrand = [&](const Real latitude) {
    const Real velocity = zonal_velocity(latitude);
    return velocity * (2.0 * rotation_rate_rad_s * std::sin(latitude) +
                       velocity * std::tan(latitude) / grid.radius_m());
  };
  std::vector<Real> balance_integral(integration_intervals + 1);
  Real previous = balance_integrand(latitude_min);
  for (std::size_t sample = 1; sample <= integration_intervals; ++sample) {
    const Real latitude = latitude_min + latitude_step * static_cast<Real>(sample);
    const Real current = balance_integrand(latitude);
    balance_integral[sample] =
        balance_integral[sample - 1] + 0.5 * latitude_step * (previous + current);
    previous = current;
  }
  const auto interpolated_balance = [&](const Real latitude) {
    const Real coordinate = std::clamp((latitude - latitude_min) / latitude_step, 0.0,
                                       static_cast<Real>(integration_intervals));
    const std::size_t lower =
        std::min(static_cast<std::size_t>(coordinate), integration_intervals - 1);
    const Real fraction = coordinate - static_cast<Real>(lower);
    return (1.0 - fraction) * balance_integral[lower] +
           fraction * balance_integral[lower + 1];
  };

  ShallowWaterState state{.time_s = time_s,
                          .step = 0,
                          .depth = std::vector<Real>(grid.cell_count()),
                          .momentum = std::vector<Vec3>(grid.cell_count())};
  Real integrated_depth = 0.0;
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
    const Vec3 position = grid.cells()[cell].center;
    const Real sine_latitude = std::clamp(dot(basis.axis, position), -1.0, 1.0);
    const Real latitude = std::asin(sine_latitude);
    const Real cosine_latitude = std::cos(latitude);
    const Real lambda = longitude(basis, position);
    const Vec3 east = cosine_latitude > 1.0e-14
                          ? cross(basis.axis, position) / cosine_latitude
                          : basis.longitude_quarter;
    Real depth = -(grid.radius_m() / gravity_m_s2) * interpolated_balance(latitude);
    if (add_height_perturbation) {
      depth +=
          perturbation_height_m * cosine_latitude *
          std::exp(-std::pow(lambda / perturbation_longitude_scale, 2)) *
          std::exp(-std::pow(
              (perturbation_latitude - latitude) / perturbation_latitude_scale, 2));
    }
    state.depth[cell] = depth;
    state.momentum[cell] = depth * zonal_velocity(latitude) * east;
    integrated_depth += grid.cells()[cell].area_m2 * depth;
  }
  const Real depth_shift = mean_depth_m - integrated_depth / grid.total_area_m2();
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
    const Vec3 velocity =
        state.depth[cell] != 0.0 ? state.momentum[cell] / state.depth[cell] : Vec3{};
    state.depth[cell] += depth_shift;
    if (!(state.depth[cell] > 0.0) || !std::isfinite(state.depth[cell])) {
      throw std::invalid_argument("Galewsky constants produce invalid depth");
    }
    state.momentum[cell] = state.depth[cell] * velocity;
  }
  return state;
}

SphericalWaveMode diagnose_depth_wave_mode(const CubedSphereGrid& grid,
                                           const ShallowWaterState& state,
                                           const Vec3 axis, const int wavenumber) {
  validate_shallow_water_state(grid, state, 0.0);
  if (wavenumber <= 0) {
    throw std::invalid_argument("wave-mode wavenumber must be positive");
  }
  const AxisBasis basis = make_axis_basis(axis);
  Real cosine_component = 0.0;
  Real sine_component = 0.0;
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
    const Real phase = wavenumber * longitude(basis, grid.cells()[cell].center);
    const Real weighted_depth = grid.cells()[cell].area_m2 * state.depth[cell];
    cosine_component += weighted_depth * std::cos(phase);
    sine_component += weighted_depth * std::sin(phase);
  }
  return {.amplitude =
              2.0 * std::hypot(cosine_component, sine_component) / grid.total_area_m2(),
          .phase_rad = std::atan2(-sine_component, cosine_component) /
                       static_cast<Real>(wavenumber)};
}

std::uint64_t shallow_water_field_checksum(const ShallowWaterState& state) noexcept {
  std::uint64_t hash = 14695981039346656037ULL;
  for (const Real depth : state.depth) {
    checksum_value(hash, depth);
  }
  for (const Vec3 momentum : state.momentum) {
    checksum_value(hash, momentum.x);
    checksum_value(hash, momentum.y);
    checksum_value(hash, momentum.z);
  }
  return hash;
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
      return make_williamson6_state(
          grid, config.run.start_time_s, config.planet.gravity_m_s2,
          config.planet.rotation_rate_rad_s, config.shallow_water.mean_depth_m,
          config.shallow_water.maximum_velocity_m_s,
          {config.shallow_water.flow_axis_x, config.shallow_water.flow_axis_y,
           config.shallow_water.flow_axis_z});
    case ShallowWaterTestCase::kGalewsky:
      return make_galewsky_state(
          grid, config.run.start_time_s, config.planet.gravity_m_s2,
          config.planet.rotation_rate_rad_s, config.shallow_water.mean_depth_m,
          config.shallow_water.maximum_velocity_m_s,
          {config.shallow_water.flow_axis_x, config.shallow_water.flow_axis_y,
           config.shallow_water.flow_axis_z});
  }
  throw std::logic_error("unknown shallow-water test case");
}

}  // namespace mps
