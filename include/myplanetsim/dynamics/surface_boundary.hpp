#pragma once

#include <filesystem>
#include <span>
#include <string>
#include <vector>

#include "myplanetsim/config/experiment_config.hpp"
#include "myplanetsim/dynamics/surface_orography.hpp"

namespace mps {

class SurfaceBoundary {
 public:
  SurfaceBoundary(std::vector<Real> surface_geopotential_m2_s2,
                  std::vector<Real> land_fraction, std::string source_id,
                  std::string source_fingerprint);

  [[nodiscard]] std::span<const Real> surface_geopotential_m2_s2() const noexcept {
    return surface_geopotential_m2_s2_;
  }
  [[nodiscard]] std::span<const Real> land_fraction() const noexcept {
    return land_fraction_;
  }
  [[nodiscard]] const std::string& source_id() const noexcept { return source_id_; }
  [[nodiscard]] const std::string& source_fingerprint() const noexcept {
    return source_fingerprint_;
  }

 private:
  std::vector<Real> surface_geopotential_m2_s2_;
  std::vector<Real> land_fraction_;
  std::string source_id_;
  std::string source_fingerprint_;
};

[[nodiscard]] SurfaceBoundary make_surface_boundary(
    const SurfaceParameters& surface, const OrographyParameters& orography,
    const CubedSphereGrid& grid, Real gravity_m_s2,
    const std::filesystem::path& source_directory = {});
[[nodiscard]] Real mixed_surface_heat_capacity(Real land_fraction,
                                               Real land_heat_capacity_j_m2_k,
                                               Real ocean_heat_capacity_j_m2_k);

}  // namespace mps
