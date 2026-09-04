#include <algorithm>
#include <cmath>
#include <filesystem>

#include "myplanetsim/diagnostics/planetary_scales.hpp"
#include "myplanetsim/dynamics/dry_hydrostatic_driver.hpp"
#include "support/test.hpp"

namespace {

[[nodiscard]] mps::ExperimentConfig earth_surface_config() {
  const auto root = std::filesystem::path(__FILE__).parent_path().parent_path();
  auto config =
      mps::load_experiment_config(root / "configs" / "phase8_earth_geography.cfg");
  config.run.end_time_s = 120.0;
  config.run.time_step_s = 30.0;
  return config;
}

}  // namespace

MPS_TEST_CASE("Earth surface fixture has bounded coast and terrain health") {
  const auto config = earth_surface_config();
  const mps::DryHydrostaticDriver driver(config);
  const auto& boundary = *driver.surface_boundary();
  bool has_ocean = false;
  bool has_land = false;
  bool has_mixed = false;
  mps::Real land_area = 0.0;
  mps::Real total_area = 0.0;
  for (std::size_t cell = 0; cell < driver.grid().cell_count(); ++cell) {
    const auto fraction = boundary.land_fraction()[cell];
    has_ocean = has_ocean || fraction == 0.0;
    has_land = has_land || fraction == 1.0;
    has_mixed = has_mixed || (fraction > 0.0 && fraction < 1.0);
    MPS_CHECK(fraction >= 0.0 && fraction <= 1.0);
    MPS_CHECK(std::isfinite(boundary.surface_geopotential_m2_s2()[cell]));
    land_area += driver.grid().cells()[cell].area_m2 * fraction;
    total_area += driver.grid().cells()[cell].area_m2;
  }
  MPS_CHECK(has_ocean);
  MPS_CHECK(has_land);
  MPS_CHECK(has_mixed);
  MPS_CHECK(land_area / total_area > 0.2);
  MPS_CHECK(land_area / total_area < 0.6);
}

MPS_TEST_CASE("surface restart is identical to uninterrupted coupling") {
  const auto config = earth_surface_config();
  const mps::DryHydrostaticDriver driver(config);
  auto uninterrupted = driver.initial_state();
  driver.advance(uninterrupted, config.run.end_time_s);

  auto restarted = driver.initial_state();
  driver.advance(restarted, 60.0);
  restarted = mps::unflatten_dry_hydrostatic_surface_state(
      restarted.time_s, restarted.step,
      mps::flatten_dry_hydrostatic_surface_state(restarted, config.vertical.levels),
      driver.grid().cell_count(), config.vertical.levels);
  driver.advance(restarted, config.run.end_time_s);
  MPS_CHECK(
      mps::flatten_dry_hydrostatic_surface_state(uninterrupted,
                                                 config.vertical.levels) ==
      mps::flatten_dry_hydrostatic_surface_state(restarted, config.vertical.levels));
}

MPS_TEST_CASE("registered dimensional similarity preserves scale metadata") {
  auto planet = mps::PlanetParameters::earth_like();
  const mps::PlanetaryScaleInputs inputs{20.0, 1.0e6, 0.01, 3.0e6, 288.0, 60.0};
  const auto base = mps::compute_planetary_scales(planet, inputs);
  constexpr mps::Real length_scale = 3.0;
  constexpr mps::Real time_scale = 2.0;
  planet.radius_m *= length_scale;
  planet.rotation_rate_rad_s /= time_scale;
  planet.gravity_m_s2 *= length_scale / (time_scale * time_scale);
  const mps::PlanetaryScaleInputs transformed{
      inputs.characteristic_wind_m_s * length_scale / time_scale,
      inputs.characteristic_length_m * length_scale,
      inputs.buoyancy_frequency_s_1 / time_scale,
      inputs.radiative_time_s * time_scale,
      inputs.reference_temperature_k * length_scale * length_scale /
          (time_scale * time_scale),
      inputs.temperature_contrast_k * length_scale * length_scale /
          (time_scale * time_scale)};
  const auto scaled = mps::compute_planetary_scales(planet, transformed);
  MPS_CHECK_NEAR(scaled.rossby_number.value, base.rossby_number.value, 1.0e-14);
  MPS_CHECK_NEAR(scaled.burger_number.value, base.burger_number.value, 1.0e-14);
  MPS_CHECK_NEAR(scaled.thermal_rossby_number.value, base.thermal_rossby_number.value,
                 1.0e-14);
  MPS_CHECK_NEAR(scaled.radiative_over_rotation.value,
                 base.radiative_over_rotation.value, 1.0e-14);
}

int main() { return mps::test::run_all(); }
