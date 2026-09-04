#include "myplanetsim/dynamics/surface_orography.hpp"

#include <cmath>
#include <numbers>
#include <stdexcept>
#include <utility>

namespace mps {

SurfaceOrography::SurfaceOrography(
    std::vector<Real> surface_geopotential_m2_s2,
    std::string source_fingerprint)
    : surface_geopotential_m2_s2_(std::move(surface_geopotential_m2_s2)),
      source_fingerprint_(std::move(source_fingerprint)) {
  if (surface_geopotential_m2_s2_.empty())
    throw std::invalid_argument("surface orography must contain cells");
  for (const auto value : surface_geopotential_m2_s2_)
    if (!std::isfinite(value))
      throw std::invalid_argument("surface orography must be finite");
  if (source_fingerprint_.empty())
    throw std::invalid_argument("surface orography fingerprint must not be empty");
}

Real dcmip_2_0_0_surface_height_m(const Vec3 position) {
  if (!is_finite(position) || !(norm_squared(position) > 0.0))
    throw std::invalid_argument("DCMIP terrain position must be finite and nonzero");
  const Vec3 unit = normalize(position);
  constexpr Vec3 mountain_centre{0.0, -1.0, 0.0};
  constexpr Real height_m = 2000.0;
  constexpr Real radius = 3.0 * std::numbers::pi_v<Real> / 4.0;
  constexpr Real half_width = std::numbers::pi_v<Real> / 16.0;
  const Real distance = safe_angle(unit, mountain_centre);
  if (!(distance < radius)) return 0.0;
  const Real envelope = 0.5 * (1.0 + std::cos(std::numbers::pi_v<Real> *
                                               distance / radius));
  const Real oscillation =
      std::cos(std::numbers::pi_v<Real> * distance / half_width);
  return height_m * envelope * oscillation * oscillation;
}

Real williamson5_surface_height_m(const Vec3 position) {
  if (!is_finite(position) || !(norm_squared(position) > 0.0))
    throw std::invalid_argument(
        "Williamson 5 terrain position must be finite and nonzero");
  const Vec3 unit = normalize(position);
  const Real longitude = std::atan2(unit.y, unit.x);
  const Real latitude = std::asin(std::clamp(unit.z, -1.0, 1.0));
  constexpr Real centre_longitude = -0.5 * std::numbers::pi_v<Real>;
  constexpr Real centre_latitude = std::numbers::pi_v<Real> / 6.0;
  constexpr Real radius = std::numbers::pi_v<Real> / 9.0;
  const Real longitude_distance =
      std::remainder(longitude - centre_longitude, 2.0 * std::numbers::pi_v<Real>);
  const Real distance =
      std::hypot(longitude_distance, latitude - centre_latitude);
  return distance < radius ? 2000.0 * (1.0 - distance / radius) : 0.0;
}

Real linear_bell_surface_height_m(const Vec3 position, const Real peak_height_m) {
  if (!is_finite(position) || !(norm_squared(position) > 0.0) ||
      !(peak_height_m > 0.0) || !std::isfinite(peak_height_m))
    throw std::invalid_argument("linear bell arguments are invalid");
  constexpr Vec3 centre{1.0, 0.0, 0.0};
  constexpr Real radius = std::numbers::pi_v<Real> / 6.0;
  const Real distance = safe_angle(normalize(position), centre);
  if (!(distance < radius)) return 0.0;
  const Real taper = std::cos(0.5 * std::numbers::pi_v<Real> * distance / radius);
  return peak_height_m * taper * taper;
}

Real jw06_surface_geopotential_m2_s2(const Vec3 position) {
  if (!is_finite(position) || !(norm_squared(position) > 0.0))
    throw std::invalid_argument("JW06 terrain position must be finite and nonzero");
  constexpr Real u0 = 35.0;
  constexpr Real eta0 = 0.252;
  constexpr Real radius_m = 6.371229e6;
  constexpr Real rotation_rate_rad_s = 7.29212e-5;
  const Real latitude = std::asin(std::clamp(normalize(position).z, -1.0, 1.0));
  const Real sine = std::sin(latitude);
  const Real cosine = std::cos(latitude);
  const Real eta_v = (1.0 - eta0) * 0.5 * std::numbers::pi_v<Real>;
  const Real vertical = std::pow(std::cos(eta_v), 1.5);
  const Real wind_shape =
      -2.0 * std::pow(sine, 6) * (cosine * cosine + 1.0 / 3.0) + 10.0 / 63.0;
  const Real rotation_shape =
      8.0 / 5.0 * std::pow(cosine, 3) * (sine * sine + 2.0 / 3.0) -
      std::numbers::pi_v<Real> / 4.0;
  return u0 * vertical *
         (u0 * vertical * wind_shape +
          radius_m * rotation_rate_rad_s * rotation_shape);
}

SurfaceOrography make_surface_orography(const OrographyParameters& parameters,
                                        const CubedSphereGrid& grid,
                                        const Real gravity_m_s2) {
  if (!(gravity_m_s2 > 0.0) || !std::isfinite(gravity_m_s2))
    throw std::invalid_argument("surface orography gravity must be positive");
  if (parameters.kind == OrographyKind::kFlat)
    return SurfaceOrography(std::vector<Real>(grid.cell_count(), 0.0),
                            "analytic:flat");
  if (parameters.kind == OrographyKind::kDcmip200) {
    std::vector<Real> geopotential(grid.cell_count());
    for (std::size_t cell = 0; cell < grid.cell_count(); ++cell)
      geopotential[cell] = gravity_m_s2 *
                           dcmip_2_0_0_surface_height_m(grid.cells()[cell].center);
    return SurfaceOrography(std::move(geopotential), "analytic:dcmip_2_0_0");
  }
  if (parameters.kind == OrographyKind::kWilliamson5) {
    std::vector<Real> geopotential(grid.cell_count());
    for (std::size_t cell = 0; cell < grid.cell_count(); ++cell)
      geopotential[cell] = gravity_m_s2 *
                           williamson5_surface_height_m(grid.cells()[cell].center);
    return SurfaceOrography(std::move(geopotential), "analytic:williamson5");
  }
  if (parameters.kind == OrographyKind::kLinearBell) {
    std::vector<Real> geopotential(grid.cell_count());
    for (std::size_t cell = 0; cell < grid.cell_count(); ++cell)
      geopotential[cell] = gravity_m_s2 *
                           linear_bell_surface_height_m(grid.cells()[cell].center);
    return SurfaceOrography(std::move(geopotential), "analytic:linear_bell");
  }
  if (parameters.kind == OrographyKind::kJw06) {
    std::vector<Real> geopotential(grid.cell_count());
    for (std::size_t cell = 0; cell < grid.cell_count(); ++cell)
      geopotential[cell] =
          jw06_surface_geopotential_m2_s2(grid.cells()[cell].center);
    return SurfaceOrography(std::move(geopotential), "analytic:jw06");
  }
  throw std::invalid_argument("selected orography initializer is not implemented");
}

}  // namespace mps
