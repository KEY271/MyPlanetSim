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

constexpr std::string_view kValidShallowWaterConfig = R"(
experiment.kind = shallow_water
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
shallow_water.test_case = williamson2
shallow_water.scheme = rusanov
shallow_water.reconstruction = linear
shallow_water.limiter = barth_jespersen
shallow_water.cfl = 0.5
shallow_water.mean_depth_m = 10
shallow_water.depth_floor_m = 0.1
shallow_water.diffusion_kind = laplacian
shallow_water.diffusion_coefficient = 0.25
shallow_water.flow_axis_x = 1
shallow_water.flow_axis_y = 2
shallow_water.flow_axis_z = 3
shallow_water.maximum_velocity_m_s = 4
diagnostics.interval_steps = 2
output.directory = output
)";

constexpr std::string_view kValidVerticalConfig = R"(
experiment.kind = vertical_column
planet.radius_m = 2
planet.rotation_rate_rad_s = 0.1
planet.gravity_m_s2 = 10
planet.gas_constant_j_kg_k = 287
planet.heat_capacity_cp_j_kg_k = 1004
planet.reference_pressure_pa = 100000
run.start_time_s = 0
run.end_time_s = 10
run.time_step_s = 1
run.random_seed = 7
vertical.test_case = dry_adiabatic
vertical.levels = 2
vertical.a_half_pa = 1000,500,0
vertical.b_half = 0,0.5,1
vertical.surface_pressure_pa = 100000
vertical.minimum_surface_pressure_pa = 90000
vertical.maximum_surface_pressure_pa = 110000
vertical.minimum_pressure_thickness_pa = 100
vertical.surface_geopotential_m2_s2 = 0
vertical.initial_temperature_k = 280
vertical.initial_potential_temperature_k = 300
vertical.temperature_floor_k = 100
vertical.transport_scheme = donor_cell
vertical.limiter = none
vertical.cfl = 0.5
vertical.forcing_amplitude = 0
diagnostics.interval_steps = 2
output.directory = output
)";

constexpr std::string_view kValidDryHydrostaticConfig = R"(
experiment.kind = dry_hydrostatic
planet.radius_m = 6371220
planet.rotation_rate_rad_s = 7.292115e-5
planet.gravity_m_s2 = 9.80616
planet.gas_constant_j_kg_k = 287
planet.heat_capacity_cp_j_kg_k = 1004
planet.reference_pressure_pa = 100000
run.start_time_s = 0
run.end_time_s = 10
run.time_step_s = 1
run.random_seed = 7
grid.cells_per_panel = 4
vertical.levels = 2
vertical.a_half_pa = 1000,500,0
vertical.b_half = 0,0.5,1
vertical.surface_pressure_pa = 100000
vertical.minimum_surface_pressure_pa = 90000
vertical.maximum_surface_pressure_pa = 110000
vertical.minimum_pressure_thickness_pa = 100
vertical.surface_geopotential_m2_s2 = 0
vertical.initial_temperature_k = 280
vertical.initial_potential_temperature_k = 300
vertical.temperature_floor_k = 100
vertical.transport_scheme = linear
vertical.limiter = minmod
vertical.cfl = 0.5
dry_hydrostatic.test_case = isothermal_rest
dry_hydrostatic.reconstruction = linear
dry_hydrostatic.limiter = barth_jespersen
dry_hydrostatic.cfl = 0.45
dry_hydrostatic.diffusion_kind = none
dry_hydrostatic.diffusion_coefficient = 0
diagnostics.interval_steps = 2
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

MPS_TEST_CASE("shallow-water configuration has a strict canonical round trip") {
  const auto first = parse(kValidShallowWaterConfig);
  MPS_CHECK(first.kind == mps::ExperimentKind::kShallowWater);
  MPS_CHECK(first.shallow_water.scheme == mps::ShallowWaterScheme::kRusanov);
  MPS_CHECK_NEAR(first.shallow_water.mean_depth_m, 10.0, 0.0);
  std::ostringstream output;
  mps::write_experiment_config(output, first);
  const auto second = parse(output.str());
  MPS_CHECK_EQ(second.shallow_water.diffusion_coefficient,
               first.shallow_water.diffusion_coefficient);
  MPS_CHECK_EQ(second.diagnostics.interval_steps, 2U);
  MPS_CHECK_THROWS_AS(
      parse(std::string(kValidShallowWaterConfig) + "transport.cfl = 0.5\n"),
      std::runtime_error);
}

MPS_TEST_CASE("shallow-water configuration rejects invalid numerical choices") {
  auto replace = [](std::string text, const std::string_view old_value,
                    const std::string_view new_value) {
    const auto position = text.find(old_value);
    text.replace(position, old_value.size(), new_value);
    return text;
  };
  MPS_CHECK_THROWS_AS(
      parse(replace(std::string(kValidShallowWaterConfig), "shallow_water.cfl = 0.5",
                    "shallow_water.cfl = 0")),
      std::invalid_argument);
  MPS_CHECK_THROWS_AS(parse(replace(std::string(kValidShallowWaterConfig),
                                    "shallow_water.depth_floor_m = 0.1",
                                    "shallow_water.depth_floor_m = 10")),
                      std::invalid_argument);
  MPS_CHECK_THROWS_AS(parse(replace(std::string(kValidShallowWaterConfig),
                                    "shallow_water.diffusion_coefficient = 0.25",
                                    "shallow_water.diffusion_coefficient = 0")),
                      std::invalid_argument);
  MPS_CHECK_THROWS_AS(parse(replace(std::string(kValidShallowWaterConfig),
                                    "diagnostics.interval_steps = 2",
                                    "diagnostics.interval_steps = 0")),
                      std::invalid_argument);
}

MPS_TEST_CASE("Williamson 5 is Rusanov-only and requires its named terrain") {
  auto text = std::string(kValidShallowWaterConfig);
  auto position = text.find("williamson2");
  text.replace(position, std::string("williamson2").size(), "williamson5");
  const auto valid = parse(text + "orography.kind = williamson5\n");
  MPS_CHECK(valid.shallow_water.test_case ==
            mps::ShallowWaterTestCase::kWilliamson5);
  position = text.find("shallow_water.scheme = rusanov");
  text.replace(position, std::string("shallow_water.scheme = rusanov").size(),
               "shallow_water.scheme = compatible");
  MPS_CHECK_THROWS_AS(parse(text + "orography.kind = williamson5\n"),
                      std::invalid_argument);
}

MPS_TEST_CASE("vertical-column configuration has a strict canonical round trip") {
  const auto first = parse(kValidVerticalConfig);
  MPS_CHECK(first.kind == mps::ExperimentKind::kVerticalColumn);
  MPS_CHECK_EQ(first.vertical.a_half_pa.size(), 3U);
  std::ostringstream output;
  mps::write_experiment_config(output, first);
  const auto second = parse(output.str());
  MPS_CHECK_NEAR(second.vertical.b_half[1], first.vertical.b_half[1], 0.0);
  MPS_CHECK_THROWS_AS(
      parse(std::string(kValidVerticalConfig) + "grid.cells_per_panel = 2\n"),
      std::runtime_error);
  MPS_CHECK_THROWS_AS(
      parse(std::string(kValidVerticalConfig) + "vertical.a_half_pa = 1,,0\n"),
      std::runtime_error);
}

MPS_TEST_CASE("dry-hydrostatic configuration reuses vertical schema strictly") {
  const auto first = parse(kValidDryHydrostaticConfig);
  MPS_CHECK(first.kind == mps::ExperimentKind::kDryHydrostatic);
  MPS_CHECK_EQ(first.grid.cells_per_panel, 4);
  MPS_CHECK(first.dry_hydrostatic.test_case ==
            mps::DryHydrostaticTestCase::kIsothermalRest);
  std::ostringstream output;
  mps::write_experiment_config(output, first);
  const auto second = parse(output.str());
  MPS_CHECK_EQ(second.vertical.levels, 2);
  MPS_CHECK_NEAR(second.dry_hydrostatic.cfl, 0.45, 0.0);
  MPS_CHECK_THROWS_AS(parse(std::string(kValidDryHydrostaticConfig) +
                            "vertical.forcing_amplitude = 0\n"),
                      std::runtime_error);
  auto dcmip_text = std::string(kValidDryHydrostaticConfig);
  const auto test_case = dcmip_text.find("isothermal_rest");
  dcmip_text.replace(test_case, std::string("isothermal_rest").size(),
                     "dcmip_2_0_0_rest");
  const auto dcmip =
      parse(dcmip_text + "orography.kind = dcmip_2_0_0\n");
  MPS_CHECK(dcmip.dry_hydrostatic.test_case ==
            mps::DryHydrostaticTestCase::kDcmip200Rest);
}

MPS_TEST_CASE("orography configuration is bounded and conditionally canonical") {
  const auto flat = parse(kValidDryHydrostaticConfig);
  std::ostringstream flat_output;
  mps::write_experiment_config(flat_output, flat);
  MPS_CHECK(flat_output.str().find("orography.") == std::string::npos);

  const auto analytic =
      parse(std::string(kValidDryHydrostaticConfig) +
            "orography.kind = linear_bell\n");
  MPS_CHECK(analytic.orography.kind == mps::OrographyKind::kLinearBell);
  std::ostringstream analytic_output;
  mps::write_experiment_config(analytic_output, analytic);
  MPS_CHECK(analytic_output.str().find("orography.kind = linear_bell\n") !=
            std::string::npos);

  const std::string imported =
      std::string(kValidDryHydrostaticConfig) +
      "orography.kind = latlon_csv\n"
      "orography.input_file = terrain.csv\n"
      "orography.input_fingerprint_fnv1a64 = 0123456789abcdef\n"
      "orography.smoothing_passes = 2\n";
  const auto csv = parse(imported);
  MPS_CHECK_EQ(csv.orography.smoothing_passes, 2);
  std::ostringstream csv_output;
  mps::write_experiment_config(csv_output, csv);
  MPS_CHECK_EQ(parse(csv_output.str()).orography.input_file, "terrain.csv");

  MPS_CHECK_THROWS_AS(
      parse(std::string(kValidDryHydrostaticConfig) +
            "orography.kind = linear_bell\n"
            "orography.smoothing_passes = 1\n"),
      std::runtime_error);
  MPS_CHECK_THROWS_AS(
      parse(std::string(kValidDryHydrostaticConfig) +
            "orography.kind = latlon_csv\n"),
      std::runtime_error);
}

int main() { return mps::test::run_all(); }
