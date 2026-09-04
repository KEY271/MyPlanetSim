#pragma once

#include <span>
#include <filesystem>
#include <string_view>
#include <string>
#include <vector>

#include "myplanetsim/config/experiment_config.hpp"
#include "myplanetsim/grid/cubed_sphere_grid.hpp"

namespace mps {

class SurfaceOrography {
 public:
  SurfaceOrography(std::vector<Real> surface_geopotential_m2_s2,
                   std::string source_fingerprint);

  [[nodiscard]] std::span<const Real> surface_geopotential_m2_s2() const noexcept {
    return surface_geopotential_m2_s2_;
  }
  [[nodiscard]] const std::string& source_fingerprint() const noexcept {
    return source_fingerprint_;
  }

 private:
  std::vector<Real> surface_geopotential_m2_s2_;
  std::string source_fingerprint_;
};

[[nodiscard]] SurfaceOrography make_surface_orography(
    const OrographyParameters& parameters, const CubedSphereGrid& grid,
    Real gravity_m_s2, const std::filesystem::path& source_directory = {});
[[nodiscard]] Real dcmip_2_0_0_surface_height_m(Vec3 position);
[[nodiscard]] Real williamson5_surface_height_m(Vec3 position);
[[nodiscard]] Real linear_bell_surface_height_m(Vec3 position,
                                                Real peak_height_m = 10.0);
[[nodiscard]] Real jw06_surface_geopotential_m2_s2(Vec3 position);
[[nodiscard]] std::string fnv1a64_hex(std::string_view bytes);
[[nodiscard]] std::vector<Real> smooth_surface_geopotential(
    const CubedSphereGrid& grid, std::span<const Real> values, Index passes);

}  // namespace mps
