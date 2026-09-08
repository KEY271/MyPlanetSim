#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "myplanetsim/physics/dry_convective_adjustment.hpp"
#include "myplanetsim/physics/simple_betts_miller.hpp"
#include "support/test.hpp"

namespace {

struct FixtureColumn {
  std::string name{};
  std::vector<double> pressure_half{};
  std::vector<double> pressure_full{};
  std::vector<double> temperature{};
  std::vector<double> vapor{};
  std::string expected_branch{};
};

[[nodiscard]] std::vector<std::string> fields(const std::string& line) {
  std::vector<std::string> result;
  std::stringstream stream(line);
  std::string field;
  while (std::getline(stream, field, ',')) result.push_back(field);
  return result;
}

[[nodiscard]] std::vector<FixtureColumn> load_fixture() {
  std::ifstream input(std::string(MPS_PHASE13_FIXTURE_DIRECTORY) +
                      "/phase13_sbm_columns.csv");
  if (!input) throw std::runtime_error("cannot open Phase 13 SBM fixture");
  std::vector<FixtureColumn> result;
  std::string line;
  while (std::getline(input, line)) {
    if (line.empty() || line.front() == '#' || line.starts_with("case,")) continue;
    const auto value = fields(line);
    if (value.size() != 8) throw std::runtime_error("invalid Phase 13 SBM fixture row");
    if (result.empty() || result.back().name != value[0]) {
      result.push_back({.name = value[0], .expected_branch = value[7]});
      result.back().pressure_half.push_back(std::stod(value[2]));
    }
    auto& column = result.back();
    column.pressure_half.push_back(std::stod(value[3]));
    column.pressure_full.push_back(std::stod(value[4]));
    column.temperature.push_back(std::stod(value[5]));
    column.vapor.push_back(std::stod(value[6]));
  }
  return result;
}

[[nodiscard]] mps::SimpleBettsMillerResult run_fixture_column(
    const FixtureColumn& column, const double time_step = 1800.0) {
  constexpr mps::DiluteMoistThermodynamics moist;
  std::vector<double> mass(column.temperature.size());
  std::vector<double> exner(column.temperature.size());
  std::vector<double> exner_half(column.pressure_half.size());
  std::vector<double> theta(column.temperature.size());
  for (std::size_t level = 0; level < mass.size(); ++level) {
    mass[level] =
        (column.pressure_half[level + 1] - column.pressure_half[level]) / 9.80665;
    exner[level] =
        std::pow(column.pressure_full[level] / 100000.0,
                 moist.gas_constant_dry_air_j_kg_k / moist.heat_capacity_cp_j_kg_k);
    theta[level] = column.temperature[level] / exner[level];
  }
  for (std::size_t level = 0; level < exner_half.size(); ++level)
    exner_half[level] =
        std::pow(column.pressure_half[level] / 100000.0,
                 moist.gas_constant_dry_air_j_kg_k / moist.heat_capacity_cp_j_kg_k);
  const auto dry = mps::dry_convective_adjustment(
      {.potential_temperature_k = theta,
       .air_mass_kg_m2 = mass,
       .exner_full = exner,
       .exner_half = exner_half,
       .heat_capacity_cp_j_kg_k = moist.heat_capacity_cp_j_kg_k,
       .heat_capacity_cv_j_kg_k =
           moist.heat_capacity_cp_j_kg_k - moist.gas_constant_dry_air_j_kg_k});
  std::vector<double> preprocessed_temperature(theta.size());
  for (std::size_t level = 0; level < theta.size(); ++level)
    preprocessed_temperature[level] =
        dry.adjusted_potential_temperature_k[level] * exner[level];
  return mps::simple_betts_miller_adjustment({.temperature_k = preprocessed_temperature,
                                              .vapor_mixing_ratio = column.vapor,
                                              .air_mass_kg_m2 = mass,
                                              .pressure_full_pa = column.pressure_full,
                                              .pressure_half_pa = column.pressure_half,
                                              .gravity_m_s2 = 9.80665,
                                              .relative_humidity_reference = 0.8,
                                              .relaxation_time_s = 7200.0,
                                              .time_step_s = time_step,
                                              .minimum_temperature_k = 150.0,
                                              .thermodynamics = moist});
}

[[nodiscard]] std::string branch_name(const mps::MoistConvectionBranch branch) {
  if (branch == mps::MoistConvectionBranch::kDeep) return "deep";
  if (branch == mps::MoistConvectionBranch::kShallow) return "shallow";
  return "none";
}

}  // namespace

MPS_TEST_CASE("fixed Phase 13 columns select none deep and shallow branches") {
  const auto fixture = load_fixture();
  MPS_CHECK_EQ(fixture.size(), 3U);
  for (const auto& column : fixture) {
    const auto result = run_fixture_column(column);
    MPS_CHECK_EQ(branch_name(result.diagnostics.branch), column.expected_branch);
  }
}

MPS_TEST_CASE("deep adjustment exports rain and conserves moist enthalpy") {
  const auto fixture = load_fixture();
  const auto result = run_fixture_column(fixture[1]);
  MPS_CHECK(result.diagnostics.branch == mps::MoistConvectionBranch::kDeep);
  MPS_CHECK(result.diagnostics.cape_j_kg > 0.0);
  MPS_CHECK(result.diagnostics.lcl_pressure_pa > fixture[1].pressure_full.front());
  MPS_CHECK(result.diagnostics.lcl_pressure_pa <= fixture[1].pressure_full.back());
  MPS_CHECK(result.diagnostics.convective_rain_kg_m2 > 0.0);
  MPS_CHECK_NEAR(result.diagnostics.column_water_change_kg_m2,
                 -result.diagnostics.convective_rain_kg_m2, 1e-13);
  MPS_CHECK_NEAR(result.diagnostics.moist_enthalpy_change_j_m2, 0.0, 1e-7);
}

MPS_TEST_CASE("shallow adjustment conserves water and exports no rain") {
  const auto fixture = load_fixture();
  const auto result = run_fixture_column(fixture[2]);
  MPS_CHECK(result.diagnostics.branch == mps::MoistConvectionBranch::kShallow);
  MPS_CHECK(result.diagnostics.lcl_pressure_pa > fixture[2].pressure_full.front());
  MPS_CHECK(result.diagnostics.lcl_pressure_pa < fixture[2].pressure_full.back());
  MPS_CHECK_EQ(result.diagnostics.convective_rain_kg_m2, 0.0);
  MPS_CHECK(result.participation_fraction[2] > 0.0);
  MPS_CHECK(result.participation_fraction[2] < 1.0);
  MPS_CHECK_NEAR(result.diagnostics.column_water_change_kg_m2, 0.0, 1e-13);
  MPS_CHECK_NEAR(result.diagnostics.moist_enthalpy_change_j_m2, 0.0, 1e-7);
}

MPS_TEST_CASE("dry neutral parcel follows the dry adiabat without false CAPE") {
  const auto fixture = load_fixture();
  const auto result = run_fixture_column(fixture[0]);
  constexpr mps::DiluteMoistThermodynamics moist;
  const double kappa =
      moist.gas_constant_dry_air_j_kg_k / moist.heat_capacity_cp_j_kg_k;
  const std::size_t bottom = result.parcel_temperature_k.size() - 1;
  for (std::size_t level = 0; level < result.parcel_temperature_k.size(); ++level) {
    const double expected =
        result.parcel_temperature_k[bottom] *
        std::pow(fixture[0].pressure_full[level] / fixture[0].pressure_full[bottom],
                 kappa);
    MPS_CHECK_NEAR(result.parcel_temperature_k[level], expected, 1e-12);
  }
  MPS_CHECK_EQ(result.diagnostics.cape_j_kg, 0.0);
}

MPS_TEST_CASE("backward Euler relaxation remains bounded for a large time step") {
  const auto fixture = load_fixture();
  const auto short_step = run_fixture_column(fixture[1], 7200.0);
  const auto long_step = run_fixture_column(fixture[1], 72000.0);
  MPS_CHECK(long_step.diagnostics.convective_rain_kg_m2 >
            short_step.diagnostics.convective_rain_kg_m2);
  for (std::size_t level = 0; level < long_step.temperature_k.size(); ++level) {
    const double target = long_step.reference_temperature_k[level];
    MPS_CHECK(std::abs(long_step.temperature_k[level] - target) <=
              std::abs(short_step.temperature_k[level] - target) + 1e-13);
    MPS_CHECK(long_step.vapor_mixing_ratio[level] >= 0.0);
  }
}

int main() { return mps::test::run_all(); }
