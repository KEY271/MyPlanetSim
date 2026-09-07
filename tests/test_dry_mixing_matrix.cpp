// Phase 12 P12.08 global comparison matrix. The four physics combinations of the plan
// share one grid and one initial state so the effect of the convective adjustment, of
// the new boundary layer, and of retiring the legacy surface exchange can be separated.
// CI runs the connection sizes (N=2/4, K=4/8); the shipped presets carry the physics
// sizes and are only started here.
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <string>
#include <vector>

#include "myplanetsim/diagnostics/dry_hydrostatic_diagnostics.hpp"
#include "myplanetsim/dynamics/dry_hydrostatic_driver.hpp"
#include "myplanetsim/vertical/hybrid_pressure_coordinate.hpp"
#include "support/test.hpp"

namespace {

constexpr mps::Real kSurfaceRefinement = 5.0;

[[nodiscard]] mps::ExperimentConfig preset(const std::string& name) {
  const auto root = std::filesystem::path(__FILE__).parent_path().parent_path();
  return mps::load_experiment_config(root / "configs" / ("phase12_" + name + ".cfg"));
}

void set_resolution(mps::ExperimentConfig& config, const int cells_per_panel,
                    const int levels) {
  config.grid.cells_per_panel = cells_per_panel;
  const auto coefficients = mps::surface_refined_sigma_coefficients(
      config.vertical.a_half_pa.front(), levels, kSurfaceRefinement);
  config.vertical.levels = levels;
  config.vertical.a_half_pa = coefficients.a_half_pa;
  config.vertical.b_half = coefficients.b_half;
}

enum class Leg { kRadiation, kConvection, kBoundaryLayer, kAll };

// The connection-sized case of the plan: a hot surface under half a day of insolation,
// so both the legacy exchange and the new boundary layer actually destabilise the
// lowest layers within a short run.
[[nodiscard]] mps::ExperimentConfig matrix_config(const Leg leg,
                                                  const int cells_per_panel = 2,
                                                  const int levels = 8) {
  auto config = preset("gray_dry_mixing");
  set_resolution(config, cells_per_panel, levels);
  config.run.end_time_s = 43200.0;
  config.run.time_step_s = 300.0;
  config.diagnostics.interval_steps = 12;
  config.surface->initial_temperature_k = 330.0;
  const bool boundary_layer = leg == Leg::kBoundaryLayer || leg == Leg::kAll;
  const bool convection = leg == Leg::kConvection || leg == Leg::kAll;
  config.boundary_layer.kind = boundary_layer ? mps::BoundaryLayerKind::kBulkKProfile
                                              : mps::BoundaryLayerKind::kNone;
  config.convection.kind =
      convection ? mps::ConvectionKind::kDryAdjustment : mps::ConvectionKind::kNone;
  // Without the boundary layer the Phase 11 owner of the sensible heat flux and of the
  // near-surface Rayleigh drag stays in place, exactly as the ownership table requires.
  config.surface->air_exchange_coefficient_w_m2_k = boundary_layer ? 0.0 : 10.0;
  config.validate();
  return config;
}

struct LegResult {
  mps::DryHydrostaticState state;
  mps::DryHydrostaticStepDiagnostics last;
  mps::Real initial_dry_mass_kg = 0.0;
  mps::Real final_dry_mass_kg = 0.0;
  mps::Real maximum_unstable_fraction_before = 0.0;
  mps::Real maximum_unstable_fraction_after = 0.0;
  mps::Real legacy_sensible_energy_j = 0.0;
  mps::Real legacy_drag_energy_j = 0.0;
  mps::Real boundary_layer_sensible_energy_j = 0.0;
  mps::Real minimum_temperature_k = 0.0;
  mps::Real maximum_temperature_k = 0.0;
};

[[nodiscard]] LegResult run(const mps::ExperimentConfig& config) {
  const mps::DryHydrostaticDriver driver(config);
  LegResult result;
  result.state = driver.initial_state();
  result.initial_dry_mass_kg =
      mps::diagnose_dry_hydrostatic_budgets(
          driver.grid(), result.state, driver.diagnose(result.state), config.planet)
          .dry_mass_kg;
  driver.advance(result.state, config.run.end_time_s,
                 [&](const mps::DryHydrostaticState&, const mps::DryHydrostaticDerived*,
                     const mps::DryHydrostaticStepDiagnostics& step) {
                   if (step.accepted_time_step_s <= 0.0) return;
                   result.last = step;
                   result.maximum_unstable_fraction_before =
                       std::max(result.maximum_unstable_fraction_before,
                                step.convection.unstable_interface_fraction_before);
                   result.maximum_unstable_fraction_after =
                       std::max(result.maximum_unstable_fraction_after,
                                step.convection.unstable_interface_fraction_after);
                   result.legacy_sensible_energy_j +=
                       step.radiation_budget.sensible_to_atmosphere_energy_j;
                   result.legacy_drag_energy_j +=
                       step.radiation_budget.rayleigh_drag_energy_j;
                   result.boundary_layer_sensible_energy_j +=
                       step.boundary_layer.sensible_to_atmosphere_energy_j;
                 });
  const auto derived = driver.diagnose(result.state);
  result.final_dry_mass_kg = mps::diagnose_dry_hydrostatic_budgets(
                                 driver.grid(), result.state, derived, config.planet)
                                 .dry_mass_kg;
  const auto [low, high] = std::ranges::minmax_element(derived.temperature_k);
  result.minimum_temperature_k = *low;
  result.maximum_temperature_k = *high;
  for (const auto value : derived.temperature_k) MPS_CHECK(std::isfinite(value));
  for (const auto value : derived.air_mass_kg_m2) MPS_CHECK(value > 0.0);
  for (const auto value : derived.tracer_mixing_ratio) MPS_CHECK(std::isfinite(value));
  return result;
}

[[nodiscard]] mps::Real difference(const mps::DryHydrostaticState& first,
                                   const mps::DryHydrostaticState& second,
                                   const int levels) {
  const auto a = mps::flatten_dry_hydrostatic_surface_state(first, levels);
  const auto b = mps::flatten_dry_hydrostatic_surface_state(second, levels);
  mps::Real result = 0.0;
  for (std::size_t index = 0; index < a.size(); ++index)
    result = std::max(
        result, std::abs(a[index] - b[index]) / std::max(1.0, std::abs(b[index])));
  return result;
}

}  // namespace

MPS_TEST_CASE(
    "the four-way matrix separates convection, boundary layer, and legacy drag") {
  const auto radiation = run(matrix_config(Leg::kRadiation));
  const auto convection = run(matrix_config(Leg::kConvection));
  const auto boundary = run(matrix_config(Leg::kBoundaryLayer));
  const auto all = run(matrix_config(Leg::kAll));

  for (const auto* leg : {&radiation, &convection, &boundary, &all}) {
    MPS_CHECK_NEAR(leg->final_dry_mass_kg, leg->initial_dry_mass_kg,
                   1e-10 * leg->initial_dry_mass_kg);
    MPS_CHECK(leg->minimum_temperature_k > 0.0);
    MPS_CHECK_EQ(leg->last.retry_count, 0U);
  }

  // Ownership: only the legs without a boundary layer keep the Phase 11 surface
  // exchange and near-surface Rayleigh drag, and only the others report boundary-layer
  // sensible heat. Nothing is counted twice.
  MPS_CHECK(radiation.legacy_sensible_energy_j != 0.0);
  MPS_CHECK(radiation.legacy_drag_energy_j != 0.0);
  MPS_CHECK_EQ(radiation.boundary_layer_sensible_energy_j, 0.0);
  MPS_CHECK(convection.legacy_sensible_energy_j != 0.0);
  MPS_CHECK(convection.legacy_drag_energy_j != 0.0);
  MPS_CHECK_EQ(convection.boundary_layer_sensible_energy_j, 0.0);
  MPS_CHECK_EQ(boundary.legacy_sensible_energy_j, 0.0);
  MPS_CHECK_EQ(boundary.legacy_drag_energy_j, 0.0);
  MPS_CHECK(boundary.boundary_layer_sensible_energy_j != 0.0);
  MPS_CHECK_EQ(all.legacy_sensible_energy_j, 0.0);
  MPS_CHECK_EQ(all.legacy_drag_energy_j, 0.0);
  MPS_CHECK(all.boundary_layer_sensible_energy_j != 0.0);

  // The instability before the adjustment is measured, not only the zero it leaves
  // behind, so the improvement is not the structural tautology of the plan.
  MPS_CHECK(convection.maximum_unstable_fraction_before > 0.0);
  MPS_CHECK_EQ(convection.maximum_unstable_fraction_after, 0.0);
  MPS_CHECK(all.maximum_unstable_fraction_before > 0.0);
  MPS_CHECK_EQ(all.maximum_unstable_fraction_after, 0.0);
  MPS_CHECK_EQ(radiation.maximum_unstable_fraction_before, 0.0);
  MPS_CHECK_EQ(boundary.maximum_unstable_fraction_before, 0.0);

  // Every process changes the trajectory, and the joint run is not either single one.
  const int levels = matrix_config(Leg::kAll).vertical.levels;
  MPS_CHECK(difference(convection.state, radiation.state, levels) > 1e-6);
  MPS_CHECK(difference(boundary.state, radiation.state, levels) > 1e-6);
  MPS_CHECK(difference(all.state, convection.state, levels) > 1e-6);
  MPS_CHECK(difference(all.state, boundary.state, levels) > 1e-6);
}

MPS_TEST_CASE("the mixing matrix connects at both CI resolutions") {
  for (const auto resolution : {2, 4}) {
    for (const auto levels : {4, 8}) {
      auto config = matrix_config(Leg::kAll, resolution, levels);
      config.run.end_time_s = 1800.0;
      config.run.time_step_s = 300.0;
      const auto leg = run(config);
      MPS_CHECK_NEAR(leg.final_dry_mass_kg, leg.initial_dry_mass_kg,
                     1e-10 * leg.initial_dry_mass_kg);
      MPS_CHECK(leg.boundary_layer_sensible_energy_j != 0.0);
      MPS_CHECK(leg.last.boundary_layer.mean_boundary_layer_height_m > 0.0);
      MPS_CHECK(std::abs(leg.last.boundary_layer.tracer_mass_change_kg) <= 1e-9);
      MPS_CHECK_EQ(leg.maximum_unstable_fraction_after, 0.0);
    }
  }
}

MPS_TEST_CASE(
    "land, ocean, and topography reach the boundary layer through the fixture") {
  auto config = preset("gray_land_ocean");
  set_resolution(config, 4, 8);
  config.run.end_time_s = 1800.0;
  config.run.time_step_s = 300.0;
  config.validate();
  const mps::DryHydrostaticDriver driver(config);
  const auto& boundary = *driver.surface_boundary();
  mps::Real minimum_land = 1.0;
  mps::Real maximum_land = 0.0;
  mps::Real maximum_geopotential = 0.0;
  for (std::size_t cell = 0; cell < driver.grid().cell_count(); ++cell) {
    minimum_land = std::min(minimum_land, boundary.land_fraction()[cell]);
    maximum_land = std::max(maximum_land, boundary.land_fraction()[cell]);
    maximum_geopotential =
        std::max(maximum_geopotential, boundary.surface_geopotential_m2_s2()[cell]);
  }
  // Fractional tiles and real topography, not a uniform surface with a flat floor.
  MPS_CHECK(minimum_land < 0.05);
  MPS_CHECK(maximum_land > 0.95);
  MPS_CHECK(maximum_geopotential > 0.0);
  const auto leg = run(config);
  MPS_CHECK_NEAR(leg.final_dry_mass_kg, leg.initial_dry_mass_kg,
                 1e-10 * leg.initial_dry_mass_kg);
  MPS_CHECK(leg.last.boundary_layer.mean_boundary_layer_height_m > 0.0);
  MPS_CHECK(std::isfinite(leg.last.boundary_layer.heat_budget_residual_j));
}

MPS_TEST_CASE(
    "the synchronous preset mixes the dayside and leaves the nightside stable") {
  auto config = preset("gray_tidally_locked");
  set_resolution(config, 4, 8);
  config.run.end_time_s = 3600.0;
  config.run.time_step_s = 300.0;
  config.surface->initial_temperature_k = 330.0;
  config.validate();
  const mps::DryHydrostaticDriver driver(config);
  auto state = driver.initial_state();
  const auto orbit =
      mps::evaluate_orbit(*config.orbit, config.planet.rotation_rate_rad_s, 0.0);
  // Cool the nightside so the two hemispheres start with opposite surface stability.
  for (std::size_t cell = 0; cell < driver.grid().cell_count(); ++cell)
    if (mps::cosine_solar_zenith(driver.grid().cells()[cell].center, orbit) == 0.0)
      state.surface_temperature_k[cell] = 250.0;
  driver.advance(state, config.run.end_time_s);
  const auto derived = driver.diagnose(state);
  const auto levels = static_cast<std::size_t>(config.vertical.levels);
  mps::Real day_bottom = 0.0;
  mps::Real night_bottom = 0.0;
  std::size_t day_cells = 0;
  std::size_t night_cells = 0;
  for (std::size_t cell = 0; cell < driver.grid().cell_count(); ++cell) {
    const mps::Real theta =
        derived.potential_temperature_k[cell * levels + levels - 1] -
        derived.potential_temperature_k[cell * levels + levels - 2];
    if (mps::cosine_solar_zenith(driver.grid().cells()[cell].center, orbit) > 0.0) {
      day_bottom += theta;
      ++day_cells;
    } else {
      night_bottom += theta;
      ++night_cells;
    }
  }
  MPS_CHECK(day_cells > 0 && night_cells > 0);
  // The heated dayside erodes the near-surface stability the cooled nightside keeps.
  MPS_CHECK(day_bottom / static_cast<mps::Real>(day_cells) >
            night_bottom / static_cast<mps::Real>(night_cells));
}

MPS_TEST_CASE(
    "coefficient sensitivity moves the surface exchange in the stated sense") {
  const auto measure = [](const mps::Real gustiness, const mps::Real critical) {
    auto config = matrix_config(Leg::kBoundaryLayer, 2, 8);
    config.boundary_layer.gustiness_m_s = gustiness;
    config.boundary_layer.critical_richardson = critical;
    config.run.end_time_s = 1800.0;
    config.validate();
    return run(config);
  };
  const auto calm = measure(0.0, 1.0);
  const auto breezy = measure(0.5, 1.0);
  const auto gusty = measure(1.0, 1.0);
  // Gustiness is the explicit unresolved-wind parameter, so raising it strengthens the
  // exchange. It is not the only source here: the run is not at rest, and the resolved
  // wind keeps `U_eff` positive even at zero gustiness. The exactly still limit is
  // gated by the zero-wind case in test_boundary_layer.cpp instead.
  MPS_CHECK(std::abs(calm.boundary_layer_sensible_energy_j) > 0.0);
  MPS_CHECK(std::abs(breezy.boundary_layer_sensible_energy_j) >
            std::abs(calm.boundary_layer_sensible_energy_j));
  MPS_CHECK(std::abs(gusty.boundary_layer_sensible_energy_j) >
            std::abs(breezy.boundary_layer_sensible_energy_j));

  const auto strict = measure(1.0, 0.25);
  const auto loose = measure(1.0, 1.0);
  MPS_CHECK(loose.last.boundary_layer.mean_boundary_layer_height_m >
            strict.last.boundary_layer.mean_boundary_layer_height_m);
}

MPS_TEST_CASE("refining the time step reduces the splitting difference") {
  const auto trajectory = [](const mps::Real dt) {
    auto config = matrix_config(Leg::kAll, 2, 8);
    config.run.end_time_s = 3600.0;
    config.run.time_step_s = dt;
    config.validate();
    return run(config).state;
  };
  const int levels = 8;
  const auto reference = trajectory(450.0 / 8.0);
  const mps::Real coarse = difference(trajectory(450.0), reference, levels);
  const mps::Real medium = difference(trajectory(225.0), reference, levels);
  const mps::Real fine = difference(trajectory(112.5), reference, levels);
  MPS_CHECK(coarse > 0.0);
  MPS_CHECK(medium < coarse);
  MPS_CHECK(fine < medium);
}

int main() { return mps::test::run_all(); }
