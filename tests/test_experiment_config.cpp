#include <filesystem>
#include <sstream>
#include <stdexcept>
#include <string>

#include "myplanetsim/config/experiment_config.hpp"
#include "myplanetsim/io/run_metadata.hpp"
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

constexpr std::string_view kValidSemiImplicitSuffix = R"(
dry_hydrostatic.time_integrator = semi_implicit
dry_hydrostatic.advective_cfl = 0.45
semi_implicit.reference_surface_pressure_pa = 100000
semi_implicit.reference_temperature_k = 280
semi_implicit.implicit_weight = 0.5
semi_implicit.wave_cfl_threshold = 0.45
semi_implicit.maximum_implicit_modes = 2
semi_implicit.nonlinear_relative_tolerance = 1e-8
semi_implicit.nonlinear_maximum_iterations = 4
semi_implicit.linear_relative_tolerance = 1e-8
semi_implicit.linear_absolute_tolerance = 1e-12
semi_implicit.linear_maximum_iterations = 40
semi_implicit.gmres_restart = 20
semi_implicit.minimum_time_step_s = 0.1
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

MPS_TEST_CASE("fractional surface configuration is optional and strict") {
  const auto old_config = parse(kValidDryHydrostaticConfig);
  MPS_CHECK(!old_config.surface.has_value());
  std::ostringstream old_text;
  mps::write_experiment_config(old_text, old_config);
  MPS_CHECK(old_text.str().find("surface.") == std::string::npos);

  const auto configured = parse(std::string(kValidDryHydrostaticConfig) +
                                "surface.geography = uniform\n"
                                "surface.uniform_land_fraction = 0.25\n");
  MPS_CHECK(configured.surface.has_value());
  MPS_CHECK_EQ(configured.surface->uniform_land_fraction, 0.25);
  std::ostringstream canonical;
  mps::write_experiment_config(canonical, configured);
  const auto round_trip = parse(canonical.str());
  MPS_CHECK_EQ(round_trip.surface->uniform_land_fraction, 0.25);

  MPS_CHECK_THROWS_AS(parse(std::string(kValidDryHydrostaticConfig) +
                            "surface.uniform_land_fraction = 0.25\n"),
                      std::runtime_error);
  MPS_CHECK_THROWS_AS(parse(std::string(kValidDryHydrostaticConfig) +
                            "surface.geography = earth\n"
                            "surface.uniform_land_fraction = 0\n"),
                      std::runtime_error);

  const auto earth =
      parse(std::string(kValidDryHydrostaticConfig) +
            "surface.geography = earth\n"
            "surface.input_file = ../data/earth/earth_surface_derived.csv\n"
            "surface.input_fingerprint_fnv1a64 = f6c9b8886a03d401\n"
            "surface.quadrature_order = 2\n"
            "surface.smoothing_passes = 0\n");
  MPS_CHECK(earth.surface->geography == mps::SurfaceGeography::kEarth);
  std::ostringstream earth_text;
  mps::write_experiment_config(earth_text, earth);
  MPS_CHECK_EQ(parse(earth_text.str()).surface->input_fingerprint_fnv1a64,
               "f6c9b8886a03d401");
}

MPS_TEST_CASE("planetary forcing has conditional orbit requirements") {
  auto base = std::string(kValidDryHydrostaticConfig);
  const auto test_case = base.find("dry_hydrostatic.test_case = isothermal_rest");
  base.replace(test_case,
               std::string("dry_hydrostatic.test_case = isothermal_rest").size(),
               "dry_hydrostatic.test_case = held_suarez");
  const auto axisymmetric = parse(base +
                                  "physics.kind = planetary_newtonian\n"
                                  "forcing.geometry = axisymmetric\n"
                                  "surface.geography = uniform\n"
                                  "surface.uniform_land_fraction = 0\n");
  MPS_CHECK(axisymmetric.physics.kind == mps::PhysicsKind::kPlanetaryNewtonian);
  MPS_CHECK(!axisymmetric.orbit.has_value());
  MPS_CHECK_THROWS_AS(parse(base + "physics.kind = planetary_newtonian\n"
                                   "forcing.geometry = substellar\n"
                                   "surface.geography = uniform\n"
                                   "surface.uniform_land_fraction = 0\n"),
                      std::runtime_error);

  const auto substellar = parse(base +
                                "physics.kind = planetary_newtonian\n"
                                "forcing.geometry = substellar\n"
                                "surface.geography = uniform\n"
                                "surface.uniform_land_fraction = 0\n"
                                "orbit.period_s = 86400\n"
                                "orbit.eccentricity = 0\n"
                                "orbit.obliquity_rad = 0\n"
                                "orbit.longitude_of_periapsis_rad = 0\n"
                                "orbit.initial_mean_anomaly_rad = 0\n"
                                "orbit.initial_substellar_longitude_rad = 0\n"
                                "star.flux_at_semimajor_axis_w_m2 = 1361\n");
  MPS_CHECK(substellar.orbit.has_value());
  std::ostringstream canonical;
  mps::write_experiment_config(canonical, substellar);
  MPS_CHECK(parse(canonical.str()).orbit.has_value());
}

MPS_TEST_CASE("registered Phase 8 presets parse and validate") {
  const auto root = std::filesystem::path(__FILE__).parent_path().parent_path();
  for (const auto* name :
       {"phase8_earth_like.cfg", "phase8_slow_rotator.cfg", "phase8_rapid_rotator.cfg",
        "phase8_tidally_locked.cfg", "phase8_earth_geography.cfg"}) {
    const auto config = mps::load_experiment_config(root / "configs" / name);
    MPS_CHECK(config.kind == mps::ExperimentKind::kDryHydrostatic);
  }
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
  MPS_CHECK(valid.shallow_water.test_case == mps::ShallowWaterTestCase::kWilliamson5);
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
  const auto dcmip = parse(dcmip_text + "orography.kind = dcmip_2_0_0\n");
  MPS_CHECK(dcmip.dry_hydrostatic.test_case ==
            mps::DryHydrostaticTestCase::kDcmip200Rest);
}

MPS_TEST_CASE("semi-implicit configuration is conditional and canonical") {
  const auto legacy = parse(kValidDryHydrostaticConfig);
  MPS_CHECK(legacy.dry_hydrostatic.time_integrator ==
            mps::DryHydrostaticTimeIntegrator::kExplicitSspRk3);
  MPS_CHECK(!legacy.semi_implicit.has_value());
  std::ostringstream legacy_text;
  mps::write_experiment_config(legacy_text, legacy);
  MPS_CHECK(legacy_text.str().find("semi_implicit") == std::string::npos);
  MPS_CHECK_EQ(mps::config_fingerprint(legacy),
               mps::config_fingerprint(parse(legacy_text.str())));

  const auto configured = parse(std::string(kValidDryHydrostaticConfig) +
                                std::string(kValidSemiImplicitSuffix));
  MPS_CHECK(configured.dry_hydrostatic.time_integrator ==
            mps::DryHydrostaticTimeIntegrator::kSemiImplicit);
  MPS_CHECK(configured.semi_implicit.has_value());
  MPS_CHECK_EQ(configured.semi_implicit->maximum_implicit_modes, 2);
  std::ostringstream canonical;
  mps::write_experiment_config(canonical, configured);
  const auto round_trip = parse(canonical.str());
  MPS_CHECK_EQ(round_trip.semi_implicit->gmres_restart, 20);
  MPS_CHECK_EQ(mps::config_fingerprint(round_trip),
               mps::config_fingerprint(configured));

  auto missing =
      std::string(kValidDryHydrostaticConfig) + std::string(kValidSemiImplicitSuffix);
  const std::string missing_line = "semi_implicit.linear_absolute_tolerance = 1e-12\n";
  missing.erase(missing.find(missing_line), missing_line.size());
  MPS_CHECK_THROWS_AS(parse(missing), std::runtime_error);

  MPS_CHECK_THROWS_AS(parse(std::string(kValidDryHydrostaticConfig) +
                            "semi_implicit.reference_temperature_k = 280\n"),
                      std::runtime_error);
  MPS_CHECK_THROWS_AS(parse(std::string(kValidShallowWaterConfig) +
                            "dry_hydrostatic.time_integrator = semi_implicit\n"),
                      std::runtime_error);
}

MPS_TEST_CASE("Held-Suarez selection is additive and strictly paired") {
  const auto none = parse(kValidDryHydrostaticConfig);
  MPS_CHECK(none.physics.kind == mps::PhysicsKind::kNone);
  std::ostringstream none_output;
  mps::write_experiment_config(none_output, none);
  MPS_CHECK(none_output.str().find("physics.") == std::string::npos);

  auto held_suarez = std::string(kValidDryHydrostaticConfig);
  auto position = held_suarez.find("planet.rotation_rate_rad_s = 7.292115e-5");
  held_suarez.replace(position,
                      std::string("planet.rotation_rate_rad_s = 7.292115e-5").size(),
                      "planet.rotation_rate_rad_s = 7.29212e-5");
  position = held_suarez.find("dry_hydrostatic.test_case = isothermal_rest");
  held_suarez.replace(position,
                      std::string("dry_hydrostatic.test_case = isothermal_rest").size(),
                      "dry_hydrostatic.test_case = held_suarez");
  held_suarez += "physics.kind = held_suarez\n";
  const auto selected = parse(held_suarez);
  MPS_CHECK(selected.physics.kind == mps::PhysicsKind::kHeldSuarez);
  MPS_CHECK(selected.dry_hydrostatic.test_case ==
            mps::DryHydrostaticTestCase::kHeldSuarez);
  std::ostringstream selected_output;
  mps::write_experiment_config(selected_output, selected);
  MPS_CHECK(selected_output.str().find("physics.kind = held_suarez\n") !=
            std::string::npos);
  MPS_CHECK(parse(selected_output.str()).physics.kind == mps::PhysicsKind::kHeldSuarez);

  MPS_CHECK_THROWS_AS(
      parse(std::string(kValidDryHydrostaticConfig) + "physics.kind = held_suarez\n"),
      std::invalid_argument);
  MPS_CHECK_THROWS_AS(parse(std::string(kValidConfig) + "physics.kind = none\n"),
                      std::runtime_error);

  position = held_suarez.find("planet.radius_m = 6371220");
  held_suarez.replace(position, std::string("planet.radius_m = 6371220").size(),
                      "planet.radius_m = 6371221");
  MPS_CHECK_THROWS_AS(parse(held_suarez), std::invalid_argument);
}

MPS_TEST_CASE("orography configuration is bounded and conditionally canonical") {
  const auto flat = parse(kValidDryHydrostaticConfig);
  std::ostringstream flat_output;
  mps::write_experiment_config(flat_output, flat);
  MPS_CHECK(flat_output.str().find("orography.") == std::string::npos);

  const auto analytic =
      parse(std::string(kValidDryHydrostaticConfig) + "orography.kind = williamson5\n");
  MPS_CHECK(analytic.orography.kind == mps::OrographyKind::kWilliamson5);
  std::ostringstream analytic_output;
  mps::write_experiment_config(analytic_output, analytic);
  MPS_CHECK(analytic_output.str().find("orography.kind = williamson5\n") !=
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

  MPS_CHECK_THROWS_AS(parse(std::string(kValidDryHydrostaticConfig) +
                            "orography.kind = linear_bell\n"
                            "orography.smoothing_passes = 1\n"),
                      std::runtime_error);
  MPS_CHECK_THROWS_AS(
      parse(std::string(kValidDryHydrostaticConfig) + "orography.kind = latlon_csv\n"),
      std::runtime_error);
  auto absolute = std::string(kValidDryHydrostaticConfig) +
                  "orography.kind = latlon_csv\n"
                  "orography.input_file = /tmp/terrain.csv\n"
                  "orography.input_fingerprint_fnv1a64 = 0123456789abcdef\n"
                  "orography.smoothing_passes = 0\n";
  MPS_CHECK_THROWS_AS(parse(absolute), std::invalid_argument);
}

int main() { return mps::test::run_all(); }
