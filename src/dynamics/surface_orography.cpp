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
  throw std::invalid_argument("selected orography initializer is not implemented");
}

}  // namespace mps
