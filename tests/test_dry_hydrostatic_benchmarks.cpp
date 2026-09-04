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
  const auto flat = mps::make_surface_orography(c.orography, grid, 10);
  for (auto kind : {mps::DryHydrostaticTestCase::kIsothermalRest,
                    mps::DryHydrostaticTestCase::kSolidBodyTransport,
                    mps::DryHydrostaticTestCase::kDcmipDeformational,
                    mps::DryHydrostaticTestCase::kDcmipHadley,
                    mps::DryHydrostaticTestCase::kLinearWave,
                    mps::DryHydrostaticTestCase::kUmjs14Steady,
                    mps::DryHydrostaticTestCase::kUmjs14Baroclinic}) {
    c.dry_hydrostatic.test_case = kind;
    auto s = mps::initialize_dry_hydrostatic_benchmark(c, grid, z, flat);
    MPS_CHECK_EQ(s.surface_pressure_pa.size(), grid.cell_count());
    for (auto v : mps::flatten_dry_hydrostatic_state(s, 2)) MPS_CHECK(std::isfinite(v));
  }
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
}
int main() { return mps::test::run_all(); }
