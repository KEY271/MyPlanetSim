#include <limits>

#include "myplanetsim/dynamics/surface_boundary.hpp"
#include "support/test.hpp"

MPS_TEST_CASE("uniform boundary adapts existing orography without thresholding") {
  const mps::CubedSphereGrid grid(3, 2.0);
  const mps::SurfaceParameters surface{.geography = mps::SurfaceGeography::kUniform,
                                       .uniform_land_fraction = 0.25};
  const auto boundary =
      mps::make_surface_boundary(surface, mps::OrographyParameters{}, grid, 10.0);
  MPS_CHECK_EQ(boundary.land_fraction().size(), grid.cell_count());
  MPS_CHECK_EQ(boundary.surface_geopotential_m2_s2().size(), grid.cell_count());
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
    MPS_CHECK_EQ(boundary.land_fraction()[cell], 0.25);
    MPS_CHECK_EQ(boundary.surface_geopotential_m2_s2()[cell], 0.0);
  }
  MPS_CHECK(boundary.source_id().starts_with("uniform+"));
  MPS_CHECK_EQ(boundary.source_fingerprint().size(), 16U);
}

MPS_TEST_CASE("mixed heat capacity is the exact fractional interpolation") {
  constexpr mps::Real land = 2.0e6;
  constexpr mps::Real ocean = 4.0e7;
  MPS_CHECK_EQ(mps::mixed_surface_heat_capacity(0.0, land, ocean), ocean);
  MPS_CHECK_EQ(mps::mixed_surface_heat_capacity(1.0, land, ocean), land);
  MPS_CHECK_EQ(mps::mixed_surface_heat_capacity(0.25, land, ocean),
               0.25 * land + 0.75 * ocean);
  MPS_CHECK_EQ(mps::mixed_surface_heat_capacity(0.5, land, ocean),
               0.5 * land + 0.5 * ocean);
  MPS_CHECK_EQ(mps::mixed_surface_heat_capacity(0.75, land, ocean),
               0.75 * land + 0.25 * ocean);
}

MPS_TEST_CASE("surface boundary rejects invalid fields") {
  MPS_CHECK_THROWS_AS(mps::SurfaceBoundary({0.0}, {0.0, 1.0}, "x", "y"),
                      std::invalid_argument);
  MPS_CHECK_THROWS_AS(
      mps::SurfaceBoundary({0.0}, {std::numeric_limits<mps::Real>::infinity()}, "x",
                           "y"),
      std::invalid_argument);
  MPS_CHECK_THROWS_AS(mps::mixed_surface_heat_capacity(1.01, 1.0, 1.0),
                      std::invalid_argument);
}

int main() { return mps::test::run_all(); }
