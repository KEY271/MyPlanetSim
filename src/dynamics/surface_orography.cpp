#include "myplanetsim/dynamics/surface_orography.hpp"

#include <cmath>
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

SurfaceOrography make_surface_orography(const OrographyParameters& parameters,
                                        const CubedSphereGrid& grid) {
  if (parameters.kind != OrographyKind::kFlat)
    throw std::invalid_argument("selected orography initializer is not implemented");
  return SurfaceOrography(std::vector<Real>(grid.cell_count(), 0.0),
                          "analytic:flat");
}

}  // namespace mps
