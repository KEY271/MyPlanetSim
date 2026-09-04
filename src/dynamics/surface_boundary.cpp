#include "myplanetsim/dynamics/surface_boundary.hpp"

#include <cmath>
#include <iomanip>
#include <limits>
#include <locale>
#include <sstream>
#include <stdexcept>
#include <utility>

#include "myplanetsim/core/validation.hpp"

namespace mps {

SurfaceBoundary::SurfaceBoundary(std::vector<Real> surface_geopotential_m2_s2,
                                 std::vector<Real> land_fraction, std::string source_id,
                                 std::string source_fingerprint)
    : surface_geopotential_m2_s2_(std::move(surface_geopotential_m2_s2)),
      land_fraction_(std::move(land_fraction)),
      source_id_(std::move(source_id)),
      source_fingerprint_(std::move(source_fingerprint)) {
  if (surface_geopotential_m2_s2_.empty() ||
      surface_geopotential_m2_s2_.size() != land_fraction_.size()) {
    throw std::invalid_argument("surface boundary fields must have equal nonzero size");
  }
  for (const Real value : surface_geopotential_m2_s2_) {
    require_finite(value, "surface geopotential");
  }
  for (const Real fraction : land_fraction_) {
    require_finite(fraction, "land fraction");
    if (fraction < 0.0 || fraction > 1.0) {
      throw std::invalid_argument("land fraction must be in [0, 1]");
    }
  }
  if (source_id_.empty() || source_fingerprint_.empty()) {
    throw std::invalid_argument("surface boundary provenance must not be empty");
  }
}

SurfaceBoundary make_surface_boundary(const SurfaceParameters& surface,
                                      const OrographyParameters& orography,
                                      const CubedSphereGrid& grid,
                                      const Real gravity_m_s2,
                                      const std::filesystem::path& source_directory) {
  require_finite(surface.uniform_land_fraction, "surface.uniform_land_fraction");
  if (surface.geography != SurfaceGeography::kUniform) {
    throw std::invalid_argument("selected surface geography is not implemented");
  }
  if (surface.uniform_land_fraction < 0.0 || surface.uniform_land_fraction > 1.0) {
    throw std::invalid_argument("surface.uniform_land_fraction must be in [0, 1]");
  }
  const auto terrain =
      make_surface_orography(orography, grid, gravity_m_s2, source_directory);
  std::ostringstream identity;
  identity.imbue(std::locale::classic());
  identity << std::setprecision(std::numeric_limits<Real>::max_digits10)
           << "uniform:" << surface.uniform_land_fraction << ':'
           << terrain.source_fingerprint();
  const std::string source_id = "uniform+" + terrain.source_fingerprint();
  return SurfaceBoundary(
      std::vector<Real>(terrain.surface_geopotential_m2_s2().begin(),
                        terrain.surface_geopotential_m2_s2().end()),
      std::vector<Real>(grid.cell_count(), surface.uniform_land_fraction), source_id,
      fnv1a64_hex(identity.str()));
}

Real mixed_surface_heat_capacity(const Real land_fraction,
                                 const Real land_heat_capacity_j_m2_k,
                                 const Real ocean_heat_capacity_j_m2_k) {
  require_finite(land_fraction, "land fraction");
  if (land_fraction < 0.0 || land_fraction > 1.0) {
    throw std::invalid_argument("land fraction must be in [0, 1]");
  }
  require_positive(land_heat_capacity_j_m2_k, "land surface heat capacity");
  require_positive(ocean_heat_capacity_j_m2_k, "ocean surface heat capacity");
  return land_fraction * land_heat_capacity_j_m2_k +
         (1.0 - land_fraction) * ocean_heat_capacity_j_m2_k;
}

}  // namespace mps
