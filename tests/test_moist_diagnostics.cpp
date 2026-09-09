#include <algorithm>
#include <cmath>
#include <vector>

#include "myplanetsim/diagnostics/moist_diagnostics.hpp"
#include "support/test.hpp"

namespace {

struct Fixture {
  mps::CubedSphereGrid grid{2, 2.0};
  mps::TracerRegistry registry{
      {{.name = "dust"}, {.name = "vapor", .role = mps::TracerRole::kWaterVapor}}};
  mps::DryHydrostaticState state;
  mps::DryHydrostaticDerived derived;
  std::vector<mps::Real> land_fraction;

  Fixture() {
    constexpr std::size_t levels = 2;
    const std::size_t cells = grid.cell_count();
    const std::size_t volume = cells * levels;
    state.tracer_count = 2;
    state.surface_pressure_pa.assign(cells, 100000.0);
    state.horizontal_momentum_mass_kg_m_s.assign(volume, {});
    state.potential_temperature_mass_k_kg_m2.assign(volume, 0.0);
    state.tracer_mass_kg_m2.assign(2 * volume, 0.0);
    state.land_water_kg_m2.resize(cells);
    land_fraction.resize(cells);
    derived.cells = cells;
    derived.levels = levels;
    derived.tracer_count = 2;
    derived.pressure_pa.resize(volume);
    derived.air_mass_kg_m2.resize(volume);
    derived.temperature_k.resize(volume);
    derived.tracer_mixing_ratio.resize(2 * volume);
    for (std::size_t cell = 0; cell < cells; ++cell) {
      land_fraction[cell] =
          static_cast<mps::Real>(cell) / static_cast<mps::Real>(cells - 1);
      state.land_water_kg_m2[cell] = land_fraction[cell] == 0.0 ? 0.0 : 10.0 + cell;
      for (std::size_t level = 0; level < levels; ++level) {
        const std::size_t scalar = mps::dry_hydrostatic_offset(cell, level, levels);
        const mps::Real mass = 1000.0 + 100.0 * level;
        const mps::Real dust = 0.001 * (1.0 + cell);
        const mps::Real vapor = 0.002 * (1.0 + level);
        derived.pressure_pa[scalar] = level == 0 ? 40000.0 : 80000.0;
        derived.air_mass_kg_m2[scalar] = mass;
        derived.temperature_k[scalar] = 270.0 + 5.0 * level;
        for (const auto& [tracer, mixing] :
             {std::pair{0U, dust}, std::pair{1U, vapor}}) {
          const std::size_t offset =
              mps::dry_hydrostatic_tracer_offset(tracer, cell, level, cells, levels);
          derived.tracer_mixing_ratio[offset] = mixing;
          state.tracer_mass_kg_m2[offset] = mass * mixing;
        }
      }
    }
  }
};

}  // namespace

MPS_TEST_CASE("tracer totals use cell area and preserve component extrema") {
  const Fixture fixture;
  const auto diagnostics = mps::diagnose_tracers(fixture.grid, fixture.state,
                                                 fixture.derived, fixture.registry);
  MPS_CHECK_EQ(diagnostics.size(), 2U);
  MPS_CHECK_NEAR(diagnostics[0].minimum_mixing_ratio, 0.001, 1e-15);
  MPS_CHECK_NEAR(diagnostics[0].maximum_mixing_ratio, 0.001 * fixture.grid.cell_count(),
                 1e-15);
  MPS_CHECK_NEAR(diagnostics[1].minimum_mixing_ratio, 0.002, 1e-15);
  MPS_CHECK_NEAR(diagnostics[1].maximum_mixing_ratio, 0.004, 1e-15);
  mps::Real direct = 0.0;
  for (std::size_t cell = 0; cell < fixture.grid.cell_count(); ++cell)
    for (std::size_t level = 0; level < 2; ++level)
      direct += fixture.grid.cells()[cell].area_m2 *
                fixture.state.tracer_mass_kg_m2[mps::dry_hydrostatic_tracer_offset(
                    1, cell, level, fixture.grid.cell_count(), 2)];
  MPS_CHECK_NEAR(diagnostics[1].mass_kg, direct, 1e-12 * direct);
}

MPS_TEST_CASE("moisture means use area-mass and area weighting") {
  const Fixture fixture;
  const mps::DiluteMoistThermodynamics thermodynamics;
  const auto diagnostics =
      mps::diagnose_moisture(fixture.grid, fixture.state, fixture.derived,
                             fixture.land_fraction, 1, thermodynamics);
  MPS_CHECK_NEAR(diagnostics.precipitable_water_kg_m2,
                 diagnostics.atmospheric_water_kg / fixture.grid.total_area_m2(),
                 1e-14);
  mps::Real expected_land_water = 0.0;
  for (std::size_t cell = 0; cell < fixture.grid.cell_count(); ++cell)
    expected_land_water += fixture.grid.cells()[cell].area_m2 *
                           fixture.land_fraction[cell] *
                           fixture.state.land_water_kg_m2[cell];
  MPS_CHECK_NEAR(diagnostics.land_water_kg, expected_land_water,
                 1e-14 * expected_land_water);
  MPS_CHECK_NEAR(diagnostics.area_mean_land_water_kg_m2,
                 expected_land_water / fixture.grid.total_area_m2(), 1e-14);
  MPS_CHECK(diagnostics.minimum_relative_humidity > 0.0);
  MPS_CHECK(diagnostics.maximum_relative_humidity >=
            diagnostics.area_mean_relative_humidity);
  MPS_CHECK_EQ(diagnostics.dilute_limit_exceedance_area_fraction, 0.0);
  MPS_CHECK_EQ(diagnostics.maximum_vapor_mixing_ratio, 0.004);
}

int main() { return mps::test::run_all(); }
