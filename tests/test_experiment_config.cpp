#include <sstream>
#include <stdexcept>
#include <string>

#include "myplanetsim/config/experiment_config.hpp"
#include "support/test.hpp"

namespace {

constexpr std::string_view kValidConfig = R"(
# comment
planet.radius_m = 6371220
planet.rotation_rate_rad_s = 7.292115e-5
planet.gravity_m_s2 = 9.80616
planet.gas_constant_j_kg_k = 287
planet.heat_capacity_cp_j_kg_k = 1004.5
planet.reference_pressure_pa = 100000
run.start_time_s = 0
run.end_time_s = 1
run.time_step_s = 0.1
run.random_seed = 42
ode.initial_value = 1
ode.decay_rate_s_1 = 2
output.directory = output data
)";

constexpr std::string_view kValidTransportConfig = R"(
experiment.kind = sphere_transport
planet.radius_m = 2
planet.rotation_rate_rad_s = 0.1
planet.gravity_m_s2 = 3
planet.gas_constant_j_kg_k = 4
planet.heat_capacity_cp_j_kg_k = 5
planet.reference_pressure_pa = 6
run.start_time_s = 0
run.end_time_s = 10
run.time_step_s = 1
run.random_seed = 7
grid.cells_per_panel = 8
grid.halo_width = 2
transport.test_case = solid_body
transport.initial_condition = gaussian_hill
transport.scheme = linear
transport.limiter = barth_jespersen
transport.cfl = 0.5
transport.rotation_axis_x = 1
transport.rotation_axis_y = 2
transport.rotation_axis_z = 3
transport.angular_speed_rad_s = 0.25
output.directory = output
)";

[[nodiscard]] mps::ExperimentConfig parse(const std::string_view text) {
  std::istringstream input{std::string(text)};
  return mps::parse_experiment_config(input);
}

}  // namespace

MPS_TEST_CASE("valid configuration is parsed") {
  const auto config = parse(kValidConfig);
  MPS_CHECK_NEAR(config.planet.radius_m, 6371220.0, 0.0);
  MPS_CHECK_NEAR(config.run.time_step_s, 0.1, 0.0);
  MPS_CHECK_EQ(config.run.random_seed, 42U);
  MPS_CHECK_NEAR(config.ode.decay_rate_s_1, 2.0, 0.0);
  MPS_CHECK_EQ(config.output_directory, "output data");
}

MPS_TEST_CASE("configuration has a canonical round trip") {
  const auto first = parse(kValidConfig);
  std::ostringstream output;
  mps::write_experiment_config(output, first);
  const auto second = parse(output.str());

  MPS_CHECK_EQ(second.planet.radius_m, first.planet.radius_m);
  MPS_CHECK_EQ(second.planet.rotation_rate_rad_s, first.planet.rotation_rate_rad_s);
  MPS_CHECK_EQ(second.run.time_step_s, first.run.time_step_s);
  MPS_CHECK_EQ(second.run.random_seed, first.run.random_seed);
  MPS_CHECK_EQ(second.ode.initial_value, first.ode.initial_value);
  MPS_CHECK_EQ(second.output_directory, first.output_directory);
}

MPS_TEST_CASE("configuration rejects unknown duplicate and missing keys") {
  MPS_CHECK_THROWS_AS(parse(std::string(kValidConfig) + "unknown.key = 1\n"),
                      std::runtime_error);
  MPS_CHECK_THROWS_AS(parse(std::string(kValidConfig) + "run.random_seed = 7\n"),
                      std::runtime_error);
  MPS_CHECK_THROWS_AS(parse("planet.radius_m = 1\n"), std::runtime_error);
}

MPS_TEST_CASE("configuration rejects invalid values") {
  auto invalid = std::string(kValidConfig);
  const auto position = invalid.find("run.time_step_s = 0.1");
  invalid.replace(position, std::string("run.time_step_s = 0.1").size(),
                  "run.time_step_s = nan");
  MPS_CHECK_THROWS_AS(parse(invalid), std::invalid_argument);

  invalid = std::string(kValidConfig);
  const auto end_position = invalid.find("run.end_time_s = 1");
  invalid.replace(end_position, std::string("run.end_time_s = 1").size(),
                  "run.end_time_s = 0");
  MPS_CHECK_THROWS_AS(parse(invalid), std::invalid_argument);
}

MPS_TEST_CASE("transport configuration has a strict canonical round trip") {
  const auto first = parse(kValidTransportConfig);
  MPS_CHECK(first.kind == mps::ExperimentKind::kSphereTransport);
  MPS_CHECK_EQ(first.grid.cells_per_panel, 8);
  MPS_CHECK(first.transport.scheme == mps::TransportScheme::kLinear);
  std::ostringstream output;
  mps::write_experiment_config(output, first);
  const auto second = parse(output.str());
  MPS_CHECK_EQ(second.transport.cfl, first.transport.cfl);
  MPS_CHECK_EQ(second.transport.rotation_axis_y, first.transport.rotation_axis_y);

  MPS_CHECK_THROWS_AS(
      parse(std::string(kValidTransportConfig) + "ode.initial_value = 1\n"),
      std::runtime_error);
}

MPS_TEST_CASE("transport configuration rejects invalid numerical choices") {
  auto invalid = std::string(kValidTransportConfig);
  auto position = invalid.find("grid.cells_per_panel = 8");
  invalid.replace(position, std::string("grid.cells_per_panel = 8").size(),
                  "grid.cells_per_panel = 0");
  MPS_CHECK_THROWS_AS(parse(invalid), std::invalid_argument);

  invalid = std::string(kValidTransportConfig);
  position = invalid.find("transport.cfl = 0.5");
  invalid.replace(position, std::string("transport.cfl = 0.5").size(),
                  "transport.cfl = 1.1");
  MPS_CHECK_THROWS_AS(parse(invalid), std::invalid_argument);

  invalid = std::string(kValidTransportConfig);
  position = invalid.find("transport.rotation_axis_x = 1");
  invalid.replace(position, std::string("transport.rotation_axis_x = 1").size(),
                  "transport.rotation_axis_x = 0");
  position = invalid.find("transport.rotation_axis_y = 2");
  invalid.replace(position, std::string("transport.rotation_axis_y = 2").size(),
                  "transport.rotation_axis_y = 0");
  position = invalid.find("transport.rotation_axis_z = 3");
  invalid.replace(position, std::string("transport.rotation_axis_z = 3").size(),
                  "transport.rotation_axis_z = 0");
  MPS_CHECK_THROWS_AS(parse(invalid), std::invalid_argument);
}

int main() { return mps::test::run_all(); }
