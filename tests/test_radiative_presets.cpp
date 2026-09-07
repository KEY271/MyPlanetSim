#include <algorithm>
#include <cmath>
#include <filesystem>
#include <string>

#include "myplanetsim/diagnostics/dry_hydrostatic_diagnostics.hpp"
#include "myplanetsim/dynamics/dry_hydrostatic_driver.hpp"
#include "support/test.hpp"

namespace {
mps::ExperimentConfig preset(const std::string& name) {
  const auto root = std::filesystem::path(__FILE__).parent_path().parent_path();
  auto config =
      mps::load_experiment_config(root / "configs" / ("phase11_" + name + ".cfg"));
  config.run.end_time_s = 120.0;
  config.run.time_step_s = 30.0;
  return config;
}
}  // namespace

MPS_TEST_CASE("radiative presets preserve mass and interface conservation") {
  for (const auto* name : {"transparent_surface", "gray_uniform", "gray_shortwave",
                           "gray_land_ocean", "gray_tidally_locked"}) {
    for (const auto resolution : {2, 4}) {
      auto config = preset(name);
      config.grid.cells_per_panel = resolution;
      const auto levels = 2 * resolution;
      config.vertical.levels = levels;
      config.vertical.a_half_pa.resize(levels + 1);
      config.vertical.b_half.resize(levels + 1);
      for (int k = 0; k <= levels; ++k) {
        const double fraction = static_cast<double>(k) / levels;
        config.vertical.a_half_pa[k] = 1000.0 * (1.0 - fraction);
        config.vertical.b_half[k] = fraction;
      }
      const mps::DryHydrostaticDriver driver(config);
      auto state = driver.initial_state();
      const auto initial = mps::diagnose_dry_hydrostatic_budgets(
          driver.grid(), state, driver.diagnose(state), config.planet);
      driver.advance(state, config.run.end_time_s);
      const auto derived = driver.diagnose(state);
      const auto final = mps::diagnose_dry_hydrostatic_budgets(driver.grid(), state,
                                                               derived, config.planet);
      MPS_CHECK_NEAR(final.dry_mass_kg, initial.dry_mass_kg,
                     1e-12 * initial.dry_mass_kg);
      for (auto temperature : derived.temperature_k)
        MPS_CHECK(std::isfinite(temperature));
      for (auto mass : derived.air_mass_kg_m2) MPS_CHECK(mass > 0.0);
      mps::DryHydrostaticRhs rhs;
      driver.rhs(state, rhs);
      MPS_CHECK(std::abs(rhs.radiation_diagnostics.interface_conservation_residual_w) <
                1e-12 * 1361.0 * driver.grid().total_area_m2());
    }
  }
}

MPS_TEST_CASE("transparent radiative trajectory matches legacy surface physics") {
  auto config = preset("transparent_surface");
  const mps::DryHydrostaticDriver radiation(config);
  auto actual = radiation.initial_state();
  radiation.advance(actual, config.run.end_time_s);
  config.physics.kind = mps::PhysicsKind::kSurfaceEnergyBalance;
  config.radiation.reset();
  const mps::DryHydrostaticDriver legacy(config);
  auto expected = legacy.initial_state();
  legacy.advance(expected, config.run.end_time_s);
  const auto first =
      mps::flatten_dry_hydrostatic_surface_state(actual, config.vertical.levels);
  const auto second =
      mps::flatten_dry_hydrostatic_surface_state(expected, config.vertical.levels);
  for (std::size_t i = 0; i < first.size(); ++i)
    MPS_CHECK_NEAR(first[i], second[i], 1e-12 * std::max(1.0, std::abs(second[i])));
}

MPS_TEST_CASE("shortwave opacity redistributes stellar absorption") {
  const mps::DryHydrostaticDriver transparent(preset("gray_uniform"));
  const mps::DryHydrostaticDriver absorbing(preset("gray_shortwave"));
  mps::DryHydrostaticRhs first, second;
  transparent.rhs(transparent.initial_state(), first);
  absorbing.rhs(absorbing.initial_state(), second);
  MPS_CHECK_NEAR(first.radiation_diagnostics.atmospheric_shortwave_heating_power_w, 0.0,
                 1.0);
  MPS_CHECK(second.radiation_diagnostics.atmospheric_shortwave_heating_power_w > 0.0);
  MPS_CHECK(second.radiation_diagnostics.surface_down_shortwave_power_w <
            first.radiation_diagnostics.surface_down_shortwave_power_w);
}

MPS_TEST_CASE("synchronous preset keeps the stellar direction fixed") {
  const auto config = preset("gray_tidally_locked");
  const auto initial =
      mps::evaluate_orbit(*config.orbit, config.planet.rotation_rate_rad_s, 0.0);
  const auto later =
      mps::evaluate_orbit(*config.orbit, config.planet.rotation_rate_rad_s, 86400.0);
  const mps::DryHydrostaticDriver driver(config);
  bool night = false;
  for (const auto& cell : driver.grid().cells()) {
    const auto cosine = mps::cosine_solar_zenith(cell.center, initial);
    MPS_CHECK_NEAR(cosine, mps::cosine_solar_zenith(cell.center, later), 1e-14);
    night = night || cosine == 0.0;
  }
  MPS_CHECK(night);
}

int main() { return mps::test::run_all(); }
