#include <algorithm>

#include "myplanetsim/dynamics/dry_hydrostatic_benchmarks.hpp"
#include "myplanetsim/dynamics/dry_hydrostatic_driver.hpp"
#include "support/test.hpp"
MPS_TEST_CASE("phase 5 benchmark presets produce shaped finite states") {
  mps::ExperimentConfig c{
      .kind = mps::ExperimentKind::kDryHydrostatic,
      .planet = {2, 0, 10, 287, 1004, 100000},
      .run = {0, 1, .1, 0},
      .grid = {2},
      .vertical = {.levels = 2,
                   .a_half_pa = {1000, 500, 0},
                   .b_half = {0, .5, 1},
                   .surface_pressure_pa = 100000,
                   .minimum_surface_pressure_pa = 90000,
                   .maximum_surface_pressure_pa = 110000,
                   .minimum_pressure_thickness_pa = 100,
                   .initial_temperature_k = 280,
                   .initial_potential_temperature_k = 300,
                   .temperature_floor_k = 100,
                   .transport_scheme = mps::VerticalTransportScheme::kDonorCell,
                   .limiter = mps::VerticalLimiterKind::kNone,
                   .cfl = .5},
      .dry_hydrostatic = {},
      .diagnostics = {1},
      .output_directory = "x"};
  mps::CubedSphereGrid grid(2, 2);
  mps::AtmosphericHybridCoordinate z({c.vertical.a_half_pa, c.vertical.b_half}, 90000,
                                     110000, 100);
  const auto flat = mps::make_surface_orography(c.orography, grid, c.planet);
  for (auto kind : {mps::DryHydrostaticTestCase::kIsothermalRest,
                    mps::DryHydrostaticTestCase::kSolidBodyTransport,
                    mps::DryHydrostaticTestCase::kDcmipDeformational,
                    mps::DryHydrostaticTestCase::kDcmipHadley,
                    mps::DryHydrostaticTestCase::kLinearWave,
                    mps::DryHydrostaticTestCase::kUmjs14Steady,
                    mps::DryHydrostaticTestCase::kUmjs14Baroclinic,
                    mps::DryHydrostaticTestCase::kHeldSuarez}) {
    c.dry_hydrostatic.test_case = kind;
    auto s = mps::initialize_dry_hydrostatic_benchmark(c, grid, z, flat);
    MPS_CHECK_EQ(s.surface_pressure_pa.size(), grid.cell_count());
    for (auto v : mps::flatten_dry_hydrostatic_state(s, 2)) MPS_CHECK(std::isfinite(v));
  }
}

MPS_TEST_CASE(
    "UMJS14 initializes a balanced thermal wind and rotational perturbation") {
  mps::ExperimentConfig config{
      .kind = mps::ExperimentKind::kDryHydrostatic,
      .planet = {6371220, 7.29212e-5, 9.80616, 287, 1004, 100000},
      .run = {0, 1, 1, 0},
      .grid = {12},
      .vertical = {.levels = 4,
                   .a_half_pa = {1000, 750, 500, 250, 0},
                   .b_half = {0, .25, .5, .75, 1},
                   .surface_pressure_pa = 100000,
                   .minimum_surface_pressure_pa = 90000,
                   .maximum_surface_pressure_pa = 110000,
                   .minimum_pressure_thickness_pa = 100,
                   .initial_temperature_k = 288,
                   .initial_potential_temperature_k = 300,
                   .temperature_floor_k = 100,
                   .transport_scheme = mps::VerticalTransportScheme::kLinear,
                   .limiter = mps::VerticalLimiterKind::kMinmod,
                   .cfl = .5},
      .dry_hydrostatic = {.test_case = mps::DryHydrostaticTestCase::kUmjs14Steady},
      .diagnostics = {1},
      .output_directory = "x"};
  const mps::CubedSphereGrid grid(config.grid.cells_per_panel, config.planet.radius_m);
  const mps::AtmosphericHybridCoordinate coordinate(
      {config.vertical.a_half_pa, config.vertical.b_half},
      config.vertical.minimum_surface_pressure_pa,
      config.vertical.maximum_surface_pressure_pa,
      config.vertical.minimum_pressure_thickness_pa);
  const auto flat = mps::make_surface_orography(config.orography, grid, config.planet);
  const auto steady =
      mps::initialize_dry_hydrostatic_benchmark(config, grid, coordinate, flat);
  config.dry_hydrostatic.test_case = mps::DryHydrostaticTestCase::kUmjs14Baroclinic;
  const auto perturbed =
      mps::initialize_dry_hydrostatic_benchmark(config, grid, coordinate, flat);

  MPS_CHECK(steady.surface_pressure_pa == perturbed.surface_pressure_pa);
  mps::Real maximum_velocity_difference = 0.0;
  for (std::size_t offset = 0; offset < steady.horizontal_momentum_mass_kg_m_s.size();
       ++offset) {
    maximum_velocity_difference =
        std::max(maximum_velocity_difference,
                 mps::norm(steady.horizontal_momentum_mass_kg_m_s[offset] -
                           perturbed.horizontal_momentum_mass_kg_m_s[offset]));
  }
  MPS_CHECK(maximum_velocity_difference > 0.0);
  MPS_CHECK(steady.potential_temperature_mass_k_kg_m2 ==
            perturbed.potential_temperature_mass_k_kg_m2);
  const auto [minimum_theta, maximum_theta] =
      std::minmax_element(steady.potential_temperature_mass_k_kg_m2.begin(),
                          steady.potential_temperature_mass_k_kg_m2.end());
  MPS_CHECK(*maximum_theta > *minimum_theta);
}

MPS_TEST_CASE("Held-Suarez initializer is seeded reproducible and level-mean free") {
  mps::ExperimentConfig config{
      .kind = mps::ExperimentKind::kDryHydrostatic,
      .planet = {6371220, 7.29212e-5, 9.80616, 287, 1004, 100000},
      .run = {0, 1, 1, 42},
      .grid = {4},
      .vertical = {.levels = 3,
                   .a_half_pa = {1000, 750, 250, 0},
                   .b_half = {0, .25, .75, 1},
                   .surface_pressure_pa = 100000,
                   .minimum_surface_pressure_pa = 90000,
                   .maximum_surface_pressure_pa = 110000,
                   .minimum_pressure_thickness_pa = 100,
                   .initial_temperature_k = 999,
                   .initial_potential_temperature_k = 300,
                   .temperature_floor_k = 100,
                   .transport_scheme = mps::VerticalTransportScheme::kLinear,
                   .limiter = mps::VerticalLimiterKind::kMinmod,
                   .cfl = .5},
      .dry_hydrostatic = {.test_case = mps::DryHydrostaticTestCase::kHeldSuarez},
      .physics = {.kind = mps::PhysicsKind::kHeldSuarez},
      .diagnostics = {1},
      .output_directory = "x"};
  config.validate();
  const mps::CubedSphereGrid grid(config.grid.cells_per_panel, config.planet.radius_m);
  const mps::AtmosphericHybridCoordinate coordinate(
      {config.vertical.a_half_pa, config.vertical.b_half},
      config.vertical.minimum_surface_pressure_pa,
      config.vertical.maximum_surface_pressure_pa,
      config.vertical.minimum_pressure_thickness_pa);
  const auto flat = mps::make_surface_orography(config.orography, grid, config.planet);
  const auto first =
      mps::initialize_dry_hydrostatic_benchmark(config, grid, coordinate, flat);
  const auto repeated =
      mps::initialize_dry_hydrostatic_benchmark(config, grid, coordinate, flat);
  MPS_CHECK(first.potential_temperature_mass_k_kg_m2 ==
            repeated.potential_temperature_mass_k_kg_m2);
  for (const auto momentum : first.horizontal_momentum_mass_kg_m_s)
    MPS_CHECK_EQ(mps::norm(momentum), 0.0);

  const auto geometry = coordinate.geometry(
      config.vertical.surface_pressure_pa, config.planet.gravity_m_s2,
      config.planet.gas_constant_j_kg_k, config.planet.heat_capacity_cp_j_kg_k,
      config.planet.reference_pressure_pa);
  for (std::size_t level = 0; level < coordinate.levels(); ++level) {
    mps::Real weighted_temperature_anomaly = 0.0;
    for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
      const auto offset = mps::dry_hydrostatic_offset(cell, level, coordinate.levels());
      const auto temperature = first.potential_temperature_mass_k_kg_m2[offset] /
                               geometry.air_mass_kg_m2[level] *
                               geometry.exner_full[level];
      MPS_CHECK(temperature > 263.9);
      MPS_CHECK(temperature < 264.1);
      weighted_temperature_anomaly +=
          grid.cells()[cell].area_m2 * (temperature - 264.0);
    }
    MPS_CHECK_NEAR(weighted_temperature_anomaly / grid.total_area_m2(), 0.0, 1e-12);
  }

  config.run.random_seed = 43;
  const auto different =
      mps::initialize_dry_hydrostatic_benchmark(config, grid, coordinate, flat);
  MPS_CHECK(first.potential_temperature_mass_k_kg_m2 !=
            different.potential_temperature_mass_k_kg_m2);
}

MPS_TEST_CASE("DCMIP 2-0-0 initializes a lapse-rate atmosphere at rest") {
  mps::ExperimentConfig config{
      .kind = mps::ExperimentKind::kDryHydrostatic,
      .planet = {6371220, 0, 9.80616, 287, 1004.5, 100000},
      .run = {0, 1, 1, 0},
      .grid = {8},
      .vertical = {.levels = 2,
                   .a_half_pa = {20544.8, 0, 0},
                   .b_half = {0, .5, 1},
                   .surface_pressure_pa = 100000,
                   .minimum_surface_pressure_pa = 75000,
                   .maximum_surface_pressure_pa = 101000,
                   .minimum_pressure_thickness_pa = 100,
                   .initial_temperature_k = 300,
                   .initial_potential_temperature_k = 300,
                   .temperature_floor_k = 100,
                   .transport_scheme = mps::VerticalTransportScheme::kDonorCell,
                   .limiter = mps::VerticalLimiterKind::kNone,
                   .cfl = .5},
      .dry_hydrostatic = {.test_case = mps::DryHydrostaticTestCase::kDcmip200Rest},
      .orography = {.kind = mps::OrographyKind::kDcmip200},
      .diagnostics = {1},
      .output_directory = "x"};
  config.validate();
  const mps::DryHydrostaticDriver driver(config);
  const auto state = driver.initial_state();
  const auto derived = driver.diagnose(state);
  const auto minimum_pressure = *std::min_element(state.surface_pressure_pa.begin(),
                                                  state.surface_pressure_pa.end());
  MPS_CHECK(minimum_pressure < 100000);
  for (const auto momentum : state.horizontal_momentum_mass_kg_m_s)
    MPS_CHECK_EQ(mps::norm(momentum), 0.0);
  for (const auto temperature : derived.temperature_k) {
    MPS_CHECK(temperature > 220);
    MPS_CHECK(temperature <= 300);
  }
  const auto rhs = driver.rhs(state);
  for (const auto momentum : rhs.tendency.momentum)
    MPS_CHECK(std::isfinite(mps::norm(momentum)));
  const auto sources = driver.diagnose_sources(derived);
  for (const auto force : sources.pressure_gradient_kg_m_s2)
    MPS_CHECK_EQ(mps::norm(force), 0.0);

  auto advanced = state;
  driver.advance(advanced, 1.0);
  const auto advanced_derived = driver.diagnose(advanced);
  for (const auto velocity : advanced_derived.velocity_m_s)
    MPS_CHECK(mps::norm(velocity) < 1.0e-12);
}
int main() { return mps::test::run_all(); }
