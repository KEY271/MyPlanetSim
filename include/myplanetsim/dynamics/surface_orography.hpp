#pragma once

#include <span>
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
    const OrographyParameters& parameters, const CubedSphereGrid& grid);

}  // namespace mps
