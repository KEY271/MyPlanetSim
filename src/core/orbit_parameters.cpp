#include "myplanetsim/core/orbit_parameters.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>

#include "myplanetsim/core/validation.hpp"

namespace mps {
namespace {

[[nodiscard]] Real wrap_signed(const Real angle) noexcept {
  constexpr Real two_pi = 2.0 * std::numbers::pi_v<Real>;
  return std::remainder(angle, two_pi);
}

[[nodiscard]] Real solve_eccentric_anomaly(const Real mean_anomaly,
                                           const Real eccentricity) {
  constexpr int maximum_iterations = 80;
  constexpr Real two_pi = 2.0 * std::numbers::pi_v<Real>;
  const Real revolutions = std::nearbyint(mean_anomaly / two_pi);
  const Real reduced_mean = mean_anomaly - revolutions * two_pi;
  Real lower = -std::numbers::pi_v<Real>;
  Real upper = std::numbers::pi_v<Real>;
  Real eccentric_anomaly = reduced_mean;
  for (int iteration = 0; iteration < maximum_iterations; ++iteration) {
    const Real residual =
        eccentric_anomaly - eccentricity * std::sin(eccentric_anomaly) - reduced_mean;
    const Real tolerance = 32.0 * std::numeric_limits<Real>::epsilon() *
                           std::max(1.0, std::abs(reduced_mean));
    if (std::abs(residual) <= tolerance) {
      return eccentric_anomaly + revolutions * two_pi;
    }
    if (residual < 0.0) {
      lower = eccentric_anomaly;
    } else {
      upper = eccentric_anomaly;
    }
    const Real derivative = 1.0 - eccentricity * std::cos(eccentric_anomaly);
    const Real candidate = eccentric_anomaly - residual / derivative;
    eccentric_anomaly =
        candidate > lower && candidate < upper ? candidate : 0.5 * (lower + upper);
  }
  throw std::runtime_error("Kepler solver failed to converge");
}

[[nodiscard]] Vec3 rotate_about_z(const Vec3 vector, const Real angle) noexcept {
  const Real cosine = std::cos(angle);
  const Real sine = std::sin(angle);
  return {cosine * vector.x - sine * vector.y, sine * vector.x + cosine * vector.y,
          vector.z};
}

}  // namespace

void OrbitParameters::validate() const {
  require_positive(period_s, "orbit.period_s");
  require_finite(eccentricity, "orbit.eccentricity");
  if (eccentricity < 0.0 || eccentricity >= 1.0) {
    throw std::invalid_argument("orbit.eccentricity must be in [0, 1)");
  }
  require_finite(obliquity_rad, "orbit.obliquity_rad");
  require_finite(longitude_of_periapsis_rad, "orbit.longitude_of_periapsis_rad");
  require_finite(initial_mean_anomaly_rad, "orbit.initial_mean_anomaly_rad");
  require_finite(initial_substellar_longitude_rad,
                 "orbit.initial_substellar_longitude_rad");
  require_positive(stellar_flux_at_semimajor_axis_w_m2,
                   "star.flux_at_semimajor_axis_w_m2");
}

Real OrbitParameters::mean_motion_rad_s() const {
  validate();
  return 2.0 * std::numbers::pi_v<Real> / period_s;
}

OrbitState evaluate_orbit(const OrbitParameters& orbit, const Real rotation_rate_rad_s,
                          const Real time_since_epoch_s) {
  orbit.validate();
  require_finite(rotation_rate_rad_s, "planet.rotation_rate_rad_s");
  require_finite(time_since_epoch_s, "time_since_epoch_s");

  const Real mean_motion = orbit.mean_motion_rad_s();
  const auto inertial_direction = [&](const Real true_anomaly) {
    const Real solar_longitude =
        true_anomaly + orbit.longitude_of_periapsis_rad + std::numbers::pi_v<Real>;
    return Vec3{std::cos(solar_longitude),
                std::cos(orbit.obliquity_rad) * std::sin(solar_longitude),
                std::sin(orbit.obliquity_rad) * std::sin(solar_longitude)};
  };
  const auto true_anomaly_at = [&](const Real time_s) {
    const Real mean_anomaly = orbit.initial_mean_anomaly_rad + mean_motion * time_s;
    const Real eccentric_anomaly =
        solve_eccentric_anomaly(mean_anomaly, orbit.eccentricity);
    return 2.0 *
           std::atan2(
               std::sqrt(1.0 + orbit.eccentricity) * std::sin(0.5 * eccentric_anomaly),
               std::sqrt(1.0 - orbit.eccentricity) * std::cos(0.5 * eccentric_anomaly));
  };

  const Vec3 initial_inertial = inertial_direction(true_anomaly_at(0.0));
  const Real initial_right_ascension =
      std::atan2(initial_inertial.y, initial_inertial.x);
  const Real rotation_phase = initial_right_ascension -
                              orbit.initial_substellar_longitude_rad +
                              rotation_rate_rad_s * time_since_epoch_s;
  const Real mean_anomaly =
      orbit.initial_mean_anomaly_rad + mean_motion * time_since_epoch_s;
  const Real eccentric_anomaly =
      solve_eccentric_anomaly(mean_anomaly, orbit.eccentricity);
  const Real true_anomaly = true_anomaly_at(time_since_epoch_s);
  const Real distance_ratio = 1.0 - orbit.eccentricity * std::cos(eccentric_anomaly);
  const Vec3 planet_fixed =
      rotate_about_z(inertial_direction(true_anomaly), -rotation_phase);
  return OrbitState{
      .star_direction_planet_fixed = normalize(planet_fixed),
      .stellar_flux_w_m2 =
          orbit.stellar_flux_at_semimajor_axis_w_m2 / (distance_ratio * distance_ratio),
      .distance_over_semimajor_axis = distance_ratio,
      .true_anomaly_rad = wrap_signed(true_anomaly),
  };
}

StellarDay stellar_day(const OrbitParameters& orbit, const Real rotation_rate_rad_s) {
  orbit.validate();
  require_finite(rotation_rate_rad_s, "planet.rotation_rate_rad_s");
  const Real difference = rotation_rate_rad_s - orbit.mean_motion_rad_s();
  const Real scale =
      std::max(std::abs(rotation_rate_rad_s), std::abs(orbit.mean_motion_rad_s()));
  const bool synchronous =
      std::abs(difference) <= 16.0 * std::numeric_limits<Real>::epsilon() * scale;
  if (synchronous) {
    return {.synchronous = true, .period_s = std::nullopt};
  }
  return {.synchronous = false,
          .period_s = 2.0 * std::numbers::pi_v<Real> / std::abs(difference)};
}

Real cosine_solar_zenith(const Vec3 unit_position, const OrbitState& state) noexcept {
  return std::max(0.0, dot(unit_position, state.star_direction_planet_fixed));
}

}  // namespace mps
