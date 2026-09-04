#include <limits>
#include <cmath>
#include <numbers>

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

MPS_TEST_CASE("Williamson 5 terrain is the published isolated cone") {
  constexpr double pi = std::numbers::pi_v<double>;
  const mps::Vec3 centre{std::cos(pi / 6) * std::cos(-pi / 2),
                         std::cos(pi / 6) * std::sin(-pi / 2),
                         std::sin(pi / 6)};
  MPS_CHECK_NEAR(mps::williamson5_surface_height_m(centre), 2000, 1e-12);
  MPS_CHECK_EQ(mps::williamson5_surface_height_m({1, 0, 0}), 0.0);
}

MPS_TEST_CASE("linear bell has fixed support and test-only amplitude scaling") {
  MPS_CHECK_NEAR(mps::linear_bell_surface_height_m({1, 0, 0}), 10, 1e-14);
  MPS_CHECK_NEAR(mps::linear_bell_surface_height_m({1, 0, 0}, 5), 5, 1e-14);
  MPS_CHECK_EQ(mps::linear_bell_surface_height_m({-1, 0, 0}), 0.0);
}

MPS_TEST_CASE("JW06 surface geopotential is smooth and zonally symmetric") {
  const double latitude = 0.7;
  const mps::Vec3 first{std::cos(latitude), 0, std::sin(latitude)};
  const mps::Vec3 second{0, std::cos(latitude), std::sin(latitude)};
  MPS_CHECK_NEAR(mps::jw06_surface_geopotential_m2_s2(first),
                 mps::jw06_surface_geopotential_m2_s2(second), 1e-11);
  MPS_CHECK(std::isfinite(mps::jw06_surface_geopotential_m2_s2({0, 0, 1})));
}

int main() { return mps::test::run_all(); }
