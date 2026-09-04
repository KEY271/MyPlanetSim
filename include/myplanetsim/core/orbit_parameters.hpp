#pragma once

#include <optional>

#include "myplanetsim/core/types.hpp"
#include "myplanetsim/geometry/vec3.hpp"

namespace mps {

struct OrbitParameters {
  Real period_s;
  Real eccentricity;
  Real obliquity_rad;
  Real longitude_of_periapsis_rad;
  Real initial_mean_anomaly_rad;
  Real initial_substellar_longitude_rad;
  Real stellar_flux_at_semimajor_axis_w_m2;

  void validate() const;
  [[nodiscard]] Real mean_motion_rad_s() const;
};

struct OrbitState {
  Vec3 star_direction_planet_fixed;
  Real stellar_flux_w_m2;
  Real distance_over_semimajor_axis;
  Real true_anomaly_rad;
};

struct StellarDay {
  bool synchronous;
  std::optional<Real> period_s;
};

[[nodiscard]] OrbitState evaluate_orbit(const OrbitParameters& orbit,
                                        Real rotation_rate_rad_s,
                                        Real time_since_epoch_s);
[[nodiscard]] StellarDay stellar_day(const OrbitParameters& orbit,
                                     Real rotation_rate_rad_s);
[[nodiscard]] Real cosine_solar_zenith(Vec3 unit_position,
                                       const OrbitState& state) noexcept;

}  // namespace mps
