#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <numbers>
#include <numeric>

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
                         std::cos(pi / 6) * std::sin(-pi / 2), std::sin(pi / 6)};
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

MPS_TEST_CASE("lat-lon CSV is fingerprinted interpolated and strict") {
  const std::string bytes =
      "longitude_deg,latitude_deg,height_m\n"
      "0,-90,10\n0,0,100\n0,90,190\n"
      "90,-90,10\n90,0,100\n90,90,190\n"
      "180,-90,10\n180,0,100\n180,90,190\n"
      "270,-90,10\n270,0,100\n270,90,190\n";
  const auto directory = std::filesystem::temp_directory_path();
  const auto path = directory / "mps_phase6_terrain.csv";
  {
    std::ofstream output(path, std::ios::binary | std::ios::trunc);
    output << bytes;
  }
  mps::OrographyParameters parameters{
      .kind = mps::OrographyKind::kLatLonCsv,
      .input_file = path.filename().string(),
      .input_fingerprint_fnv1a64 = mps::fnv1a64_hex(bytes),
      .smoothing_passes = 0};
  const mps::CubedSphereGrid grid(4, 2);
  const auto terrain = mps::make_surface_orography(parameters, grid, 10, directory);
  for (const auto value : terrain.surface_geopotential_m2_s2()) {
    MPS_CHECK(value >= 100);
    MPS_CHECK(value <= 1900);
  }
  parameters.input_fingerprint_fnv1a64 = "0000000000000000";
  MPS_CHECK_THROWS_AS(mps::make_surface_orography(parameters, grid, 10, directory),
                      std::runtime_error);

  const auto reject = [&](const std::string& invalid) {
    {
      std::ofstream output(path, std::ios::binary | std::ios::trunc);
      output << invalid;
    }
    parameters.input_fingerprint_fnv1a64 = mps::fnv1a64_hex(invalid);
    MPS_CHECK_THROWS_AS(mps::make_surface_orography(parameters, grid, 10, directory),
                        std::runtime_error);
  };
  reject(bytes + "270,90,190\n");
  reject(
      "longitude_deg,latitude_deg,height_m\n"
      "0,-90,10\n0,90,190\n90,-90,10\n");
  reject(
      "longitude_deg,latitude_deg,height_m\n"
      "0,-90,nan\n0,90,190\n90,-90,10\n90,90,190\n");
  std::filesystem::remove(path);
}

MPS_TEST_CASE("unique-edge smoothing preserves constants mean and extrema") {
  const mps::CubedSphereGrid grid(5, 2);
  const std::vector<mps::Real> constant(grid.cell_count(), 7);
  const auto same = mps::smooth_surface_geopotential(grid, constant, 3);
  MPS_CHECK(same == constant);
  std::vector<mps::Real> field(grid.cell_count());
  field.front() = 10;
  const auto smooth = mps::smooth_surface_geopotential(grid, field, 2);
  double before = 0;
  double after = 0;
  for (std::size_t cell = 0; cell < field.size(); ++cell) {
    before += grid.cells()[cell].area_m2 * field[cell];
    after += grid.cells()[cell].area_m2 * smooth[cell];
    MPS_CHECK(smooth[cell] >= 0);
    MPS_CHECK(smooth[cell] <= 10);
  }
  MPS_CHECK_NEAR(after, before, 2e-14 * before);
}

int main() { return mps::test::run_all(); }
