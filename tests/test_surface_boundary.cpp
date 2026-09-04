#include <filesystem>
#include <fstream>
#include <limits>
#include <sstream>

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

MPS_TEST_CASE("Earth product is fingerprinted and retains mixed fractions") {
  std::ostringstream bytes;
  bytes << "longitude_deg,latitude_deg,mean_surface_height_m,land_fraction\n";
  for (const double longitude : {-135.0, -45.0, 45.0, 135.0}) {
    for (const double latitude : {-90.0, 0.0, 90.0}) {
      const double fraction = longitude < 0.0 ? 0.0 : 1.0;
      bytes << longitude << ',' << latitude << ',' << 100.0 * fraction << ','
            << fraction << '\n';
    }
  }
  const auto directory = std::filesystem::temp_directory_path();
  const auto path = directory / "mps_phase8_earth_surface.csv";
  {
    std::ofstream output(path, std::ios::binary | std::ios::trunc);
    output << bytes.str();
  }
  const mps::SurfaceParameters surface{
      .geography = mps::SurfaceGeography::kEarth,
      .input_file = path.filename().string(),
      .input_fingerprint_fnv1a64 = mps::fnv1a64_hex(bytes.str()),
      .quadrature_order = 2,
      .smoothing_passes = 0};
  const mps::CubedSphereGrid grid(4, 2.0);
  const auto boundary = mps::make_surface_boundary(surface, mps::OrographyParameters{},
                                                   grid, 10.0, directory);
  bool has_ocean = false;
  bool has_land = false;
  bool has_mixed = false;
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
    const auto fraction = boundary.land_fraction()[cell];
    has_ocean = has_ocean || fraction == 0.0;
    has_land = has_land || fraction == 1.0;
    has_mixed = has_mixed || (fraction > 0.0 && fraction < 1.0);
    if (fraction == 0.0) MPS_CHECK_EQ(boundary.surface_geopotential_m2_s2()[cell], 0.0);
  }
  MPS_CHECK(has_ocean);
  MPS_CHECK(has_land);
  MPS_CHECK(has_mixed);
  MPS_CHECK_EQ(boundary.source_fingerprint(), mps::fnv1a64_hex(bytes.str()));
  std::filesystem::remove(path);
}

int main() { return mps::test::run_all(); }
