#include <limits>

#include "myplanetsim/dynamics/surface_orography.hpp"
#include "support/test.hpp"

MPS_TEST_CASE("flat surface orography owns one immutable cell field") {
  const mps::CubedSphereGrid grid(2, 3);
  const auto orography =
      mps::make_surface_orography(mps::OrographyParameters{}, grid, 10);
  MPS_CHECK_EQ(orography.surface_geopotential_m2_s2().size(), grid.cell_count());
  for (const auto value : orography.surface_geopotential_m2_s2())
    MPS_CHECK_EQ(value, 0.0);
  MPS_CHECK_EQ(orography.source_fingerprint(), "analytic:flat");
  MPS_CHECK_THROWS_AS(
      mps::SurfaceOrography({0, std::numeric_limits<mps::Real>::infinity()}, "x"),
      std::invalid_argument);
}

MPS_TEST_CASE("DCMIP 2-0-0 terrain follows the published Schar profile") {
  MPS_CHECK_NEAR(mps::dcmip_2_0_0_surface_height_m({0, -1, 0}), 2000, 1e-12);
  MPS_CHECK_EQ(mps::dcmip_2_0_0_surface_height_m({0, 1, 0}), 0.0);
  mps::OrographyParameters parameters;
  parameters.kind = mps::OrographyKind::kDcmip200;
  const mps::CubedSphereGrid grid(8, 6371220);
  const auto orography = mps::make_surface_orography(parameters, grid, 9.80616);
  for (const auto geopotential : orography.surface_geopotential_m2_s2()) {
    MPS_CHECK(geopotential >= 0.0);
    MPS_CHECK(geopotential <= 2000.0 * 9.80616);
  }
}

int main() { return mps::test::run_all(); }
