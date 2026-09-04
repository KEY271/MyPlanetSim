#include <limits>

#include "myplanetsim/dynamics/surface_orography.hpp"
#include "support/test.hpp"

MPS_TEST_CASE("flat surface orography owns one immutable cell field") {
  const mps::CubedSphereGrid grid(2, 3);
  const auto orography =
      mps::make_surface_orography(mps::OrographyParameters{}, grid);
  MPS_CHECK_EQ(orography.surface_geopotential_m2_s2().size(), grid.cell_count());
  for (const auto value : orography.surface_geopotential_m2_s2())
    MPS_CHECK_EQ(value, 0.0);
  MPS_CHECK_EQ(orography.source_fingerprint(), "analytic:flat");
  MPS_CHECK_THROWS_AS(
      mps::SurfaceOrography({0, std::numeric_limits<mps::Real>::infinity()}, "x"),
      std::invalid_argument);
}

int main() { return mps::test::run_all(); }
