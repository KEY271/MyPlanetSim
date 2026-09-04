#include <algorithm>
#include <cmath>
#include <sstream>

#include "myplanetsim/diagnostics/climate_statistics.hpp"
#include "support/test.hpp"

namespace {

struct Fixture {
  mps::CubedSphereGrid grid{2, 2.0};
  mps::DryHydrostaticState state;
  mps::DryHydrostaticDerived derived;

  Fixture() {
    const std::size_t cells = grid.cell_count();
    constexpr std::size_t levels = 2;
    state.surface_pressure_pa.assign(cells, 100000.0);
    state.horizontal_momentum_mass_kg_m_s.assign(cells * levels, {});
    state.potential_temperature_mass_k_kg_m2.assign(cells * levels, 0.0);
    state.tracer_mass_kg_m2.assign(cells * levels, 0.0);
    derived.cells = cells;
    derived.levels = levels;
    derived.pressure_pa.resize(cells * levels);
    derived.air_mass_kg_m2.resize(cells * levels);
    derived.velocity_m_s.resize(cells * levels);
    derived.temperature_k.resize(cells * levels);
    for (std::size_t cell = 0; cell < cells; ++cell) {
      const auto position = grid.cells()[cell].center;
      const mps::Vec3 east = mps::normalize(mps::Vec3{-position.y, position.x, 0.0});
      const mps::Vec3 north = mps::cross(position, east);
      for (std::size_t level = 0; level < levels; ++level) {
        const auto offset = mps::dry_hydrostatic_offset(cell, level, levels);
        derived.pressure_pa[offset] = level == 0 ? 25000.0 : 75000.0;
        derived.air_mass_kg_m2[offset] = 1000.0 + 100.0 * level;
        derived.temperature_k[offset] = 250.0 + 10.0 * level;
        derived.velocity_m_s[offset] = 3.0 * east + 2.0 * north;
      }
    }
  }
};

}  // namespace

MPS_TEST_CASE("constant fields have exact means and zero eddy moments") {
  Fixture fixture;
  mps::ClimateStatisticsAccumulator statistics(fixture.grid, 2, 0.0);
  statistics.observe(fixture.state, fixture.derived);
  fixture.state.time_s = 10.0;
  statistics.observe(fixture.state, fixture.derived);
  const auto rows = statistics.rows();
  MPS_CHECK_EQ(rows.size(), 8U);
  for (const auto& row : rows) {
    MPS_CHECK_NEAR(row.mean_pressure_pa, row.level == 0 ? 25000.0 : 75000.0, 1e-10);
    MPS_CHECK_NEAR(row.mean_temperature_k, row.level == 0 ? 250.0 : 260.0, 1e-12);
    MPS_CHECK_NEAR(row.mean_zonal_wind_m_s, 3.0, 1e-12);
    MPS_CHECK_NEAR(row.mean_meridional_wind_m_s, 2.0, 1e-12);
    MPS_CHECK_NEAR(row.eddy_momentum_flux_m2_s2, 0.0, 1e-12);
    MPS_CHECK_NEAR(row.eddy_heat_flux_k_m_s, 0.0, 1e-12);
    MPS_CHECK_NEAR(row.eddy_kinetic_energy_m2_s2, 0.0, 1e-12);
    MPS_CHECK_NEAR(row.temperature_variance_k2, 0.0, 1e-10);
    MPS_CHECK_EQ(row.accumulated_time_s, 10.0);
  }
}

MPS_TEST_CASE("online accumulation is invariant to constant-field chunking") {
  Fixture fixture;
  mps::ClimateStatisticsAccumulator one_chunk(fixture.grid, 2, 0.0);
  one_chunk.observe(fixture.state, fixture.derived);
  fixture.state.time_s = 12.0;
  one_chunk.observe(fixture.state, fixture.derived);

  fixture.state.time_s = 0.0;
  mps::ClimateStatisticsAccumulator three_chunks(fixture.grid, 2, 0.0);
  three_chunks.observe(fixture.state, fixture.derived);
  for (const auto time : {3.0, 7.0, 12.0}) {
    fixture.state.time_s = time;
    three_chunks.observe(fixture.state, fixture.derived);
  }
  const auto first = one_chunk.rows();
  const auto second = three_chunks.rows();
  MPS_CHECK_EQ(first.size(), second.size());
  for (std::size_t index = 0; index < first.size(); ++index) {
    MPS_CHECK_NEAR(first[index].mean_temperature_k, second[index].mean_temperature_k,
                   1e-13);
    MPS_CHECK_NEAR(first[index].mass_streamfunction_kg_s,
                   second[index].mass_streamfunction_kg_s, 1e-10);
    MPS_CHECK_EQ(first[index].accumulated_time_s, second[index].accumulated_time_s);
  }
}

MPS_TEST_CASE("statistics respect the spin-up boundary and deterministic CSV order") {
  Fixture fixture;
  mps::ClimateStatisticsAccumulator statistics(fixture.grid, 2, 5.0);
  statistics.observe(fixture.state, fixture.derived);
  fixture.state.time_s = 4.0;
  statistics.observe(fixture.state, fixture.derived);
  fixture.state.time_s = 8.0;
  statistics.observe(fixture.state, fixture.derived);
  const auto rows = statistics.rows();
  for (const auto& row : rows) MPS_CHECK_EQ(row.accumulated_time_s, 3.0);
  for (std::size_t index = 1; index < rows.size(); ++index) {
    MPS_CHECK(rows[index - 1].latitude_lower_deg <= rows[index].latitude_lower_deg);
    if (rows[index - 1].latitude_lower_deg == rows[index].latitude_lower_deg)
      MPS_CHECK_EQ(rows[index].level, rows[index - 1].level + 1);
  }
  std::ostringstream output;
  mps::write_climate_statistics_csv(output, rows);
  MPS_CHECK(output.str().starts_with(
      "latitude_lower_deg,latitude_upper_deg,level,mean_pressure_pa,"));
  const auto csv = output.str();
  MPS_CHECK_EQ(static_cast<std::size_t>(std::count(csv.begin(), csv.end(), '\n')),
               rows.size() + 1);
}

int main() { return mps::test::run_all(); }
