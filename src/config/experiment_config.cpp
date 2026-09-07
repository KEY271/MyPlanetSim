#include "myplanetsim/config/experiment_config.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <charconv>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <locale>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>

#include "myplanetsim/core/validation.hpp"
#include "myplanetsim/dynamics/jw06_parameters.hpp"

namespace mps {
namespace {

constexpr std::array<std::string_view, 11> kCommonRequiredKeys{
    "planet.radius_m",
    "planet.rotation_rate_rad_s",
    "planet.gravity_m_s2",
    "planet.gas_constant_j_kg_k",
    "planet.heat_capacity_cp_j_kg_k",
    "planet.reference_pressure_pa",
    "run.start_time_s",
    "run.end_time_s",
    "run.time_step_s",
    "run.random_seed",
    "output.directory",
};

constexpr std::array<std::string_view, 7> kOrbitRequiredKeys{
    "orbit.period_s",
    "orbit.eccentricity",
    "orbit.obliquity_rad",
    "orbit.longitude_of_periapsis_rad",
    "orbit.initial_mean_anomaly_rad",
    "orbit.initial_substellar_longitude_rad",
    "star.flux_at_semimajor_axis_w_m2",
};

constexpr std::array<std::string_view, 7> kSurfacePhysicsRequiredKeys{
    "surface.land_heat_capacity_j_m2_k",
    "surface.ocean_heat_capacity_j_m2_k",
    "surface.initial_temperature_k",
    "surface.albedo",
    "surface.emissivity",
    "surface.air_exchange_coefficient_w_m2_k",
    "surface.internal_heat_flux_w_m2",
};

constexpr std::array<std::string_view, 7> kRadiationRequiredKeys{
    "radiation.shortwave_absorption_m2_kg",
    "radiation.longwave_absorption_ref_m2_kg",
    "radiation.reference_pressure_pa",
    "radiation.longwave_pressure_exponent",
    "radiation.longwave_diffusivity_factor",
    "radiation.shortwave_diffuse_factor",
    "radiation.cfl",
};

constexpr std::array<std::string_view, 2> kOdeRequiredKeys{"ode.initial_value",
                                                           "ode.decay_rate_s_1"};

constexpr std::array<std::string_view, 11> kTransportRequiredKeys{
    "experiment.kind",
    "grid.cells_per_panel",
    "transport.test_case",
    "transport.initial_condition",
    "transport.scheme",
    "transport.limiter",
    "transport.cfl",
    "transport.rotation_axis_x",
    "transport.rotation_axis_y",
    "transport.rotation_axis_z",
    "transport.angular_speed_rad_s",
};

constexpr std::array<std::string_view, 16> kShallowWaterRequiredKeys{
    "experiment.kind",
    "grid.cells_per_panel",
    "shallow_water.test_case",
    "shallow_water.scheme",
    "shallow_water.reconstruction",
    "shallow_water.limiter",
    "shallow_water.cfl",
    "shallow_water.mean_depth_m",
    "shallow_water.depth_floor_m",
    "shallow_water.diffusion_kind",
    "shallow_water.diffusion_coefficient",
    "shallow_water.flow_axis_x",
    "shallow_water.flow_axis_y",
    "shallow_water.flow_axis_z",
    "shallow_water.maximum_velocity_m_s",
    "diagnostics.interval_steps",
};

constexpr std::array<std::string_view, 18> kVerticalRequiredKeys{
    "experiment.kind",
    "vertical.test_case",
    "vertical.levels",
    "vertical.a_half_pa",
    "vertical.b_half",
    "vertical.surface_pressure_pa",
    "vertical.minimum_surface_pressure_pa",
    "vertical.maximum_surface_pressure_pa",
    "vertical.minimum_pressure_thickness_pa",
    "vertical.surface_geopotential_m2_s2",
    "vertical.initial_temperature_k",
    "vertical.initial_potential_temperature_k",
    "vertical.temperature_floor_k",
    "vertical.transport_scheme",
    "vertical.limiter",
    "vertical.cfl",
    "vertical.forcing_amplitude",
    "diagnostics.interval_steps",
};

constexpr std::array<std::string_view, 24> kDryHydrostaticRequiredKeys{
    "experiment.kind",
    "grid.cells_per_panel",
    "vertical.levels",
    "vertical.a_half_pa",
    "vertical.b_half",
    "vertical.surface_pressure_pa",
    "vertical.minimum_surface_pressure_pa",
    "vertical.maximum_surface_pressure_pa",
    "vertical.minimum_pressure_thickness_pa",
    "vertical.surface_geopotential_m2_s2",
    "vertical.initial_temperature_k",
    "vertical.initial_potential_temperature_k",
    "vertical.temperature_floor_k",
    "vertical.transport_scheme",
    "vertical.limiter",
    "vertical.cfl",
    "dry_hydrostatic.test_case",
    "dry_hydrostatic.reconstruction",
    "dry_hydrostatic.limiter",
    "dry_hydrostatic.cfl",
    "dry_hydrostatic.diffusion_kind",
    "dry_hydrostatic.diffusion_coefficient",
    "diagnostics.interval_steps",
    "output.directory"};

constexpr std::array<std::string_view, 12> kSemiImplicitRequiredKeys{
    "semi_implicit.reference_surface_pressure_pa",
    "semi_implicit.reference_temperature_k",
    "semi_implicit.reference_update",
    "semi_implicit.implicit_weight",
    "semi_implicit.wave_cfl_threshold",
    "semi_implicit.maximum_implicit_modes",
    "semi_implicit.nonlinear_iterations",
    "semi_implicit.linear_relative_tolerance",
    "semi_implicit.linear_absolute_tolerance",
    "semi_implicit.linear_maximum_iterations",
    "semi_implicit.gmres_restart",
    "semi_implicit.minimum_time_step_s"};

[[nodiscard]] std::string_view trim(const std::string_view value) {
  std::size_t first = 0;
  while (first < value.size() &&
         std::isspace(static_cast<unsigned char>(value[first])) != 0) {
    ++first;
  }

  std::size_t last = value.size();
  while (last > first &&
         std::isspace(static_cast<unsigned char>(value[last - 1])) != 0) {
    --last;
  }
  return value.substr(first, last - first);
}

[[nodiscard]] std::runtime_error parse_error(const std::size_t line,
                                             const std::string_view message) {
  return std::runtime_error("configuration line " + std::to_string(line) + ": " +
                            std::string(message));
}

[[nodiscard]] Real parse_real(const std::string_view text, const std::size_t line,
                              const std::string_view key) {
  Real value = 0.0;
  const auto result = std::from_chars(text.data(), text.data() + text.size(), value);
  if (result.ec != std::errc{} || result.ptr != text.data() + text.size()) {
    throw parse_error(line, "invalid floating-point value for " + std::string(key));
  }
  return value;
}

[[nodiscard]] Seed parse_seed(const std::string_view text, const std::size_t line,
                              const std::string_view key) {
  Seed value = 0;
  const auto result = std::from_chars(text.data(), text.data() + text.size(), value);
  if (result.ec != std::errc{} || result.ptr != text.data() + text.size()) {
    throw parse_error(line, "invalid unsigned integer value for " + std::string(key));
  }
  return value;
}

[[nodiscard]] Index parse_index(const std::string_view text, const std::size_t line,
                                const std::string_view key) {
  Index value = 0;
  const auto result = std::from_chars(text.data(), text.data() + text.size(), value);
  if (result.ec != std::errc{} || result.ptr != text.data() + text.size()) {
    throw parse_error(line, "invalid integer value for " + std::string(key));
  }
  return value;
}

[[nodiscard]] std::vector<Real> parse_real_list(const std::string_view text,
                                                const std::size_t line,
                                                const std::string_view key) {
  std::vector<Real> values;
  std::size_t begin = 0;
  while (begin < text.size()) {
    const auto comma = text.find(',', begin);
    const auto item =
        trim(text.substr(begin, comma == std::string_view::npos ? std::string_view::npos
                                                                : comma - begin));
    if (item.empty()) {
      throw parse_error(line, "empty list element for " + std::string(key));
    }
    const auto value = parse_real(item, line, key);
    if (!std::isfinite(value)) {
      throw parse_error(line, "non-finite list element for " + std::string(key));
    }
    values.push_back(value);
    if (comma == std::string_view::npos) {
      break;
    }
    begin = comma + 1;
    if (begin == text.size()) {
      throw parse_error(line, "trailing comma for " + std::string(key));
    }
  }
  return values;
}

void assign_value(ExperimentConfig& config, const std::string_view key,
                  const std::string_view value, const std::size_t line) {
  if (key == "experiment.kind") {
    if (value == "ode") {
      config.kind = ExperimentKind::kOde;
    } else if (value == "sphere_transport") {
      config.kind = ExperimentKind::kSphereTransport;
    } else if (value == "shallow_water") {
      config.kind = ExperimentKind::kShallowWater;
    } else if (value == "vertical_column") {
      config.kind = ExperimentKind::kVerticalColumn;
    } else if (value == "dry_hydrostatic") {
      config.kind = ExperimentKind::kDryHydrostatic;
    } else {
      throw parse_error(line, "unknown experiment.kind " + std::string(value));
    }
  } else if (key == "planet.radius_m") {
    config.planet.radius_m = parse_real(value, line, key);
  } else if (key == "planet.rotation_rate_rad_s") {
    config.planet.rotation_rate_rad_s = parse_real(value, line, key);
  } else if (key == "planet.gravity_m_s2") {
    config.planet.gravity_m_s2 = parse_real(value, line, key);
  } else if (key == "planet.gas_constant_j_kg_k") {
    config.planet.gas_constant_j_kg_k = parse_real(value, line, key);
  } else if (key == "planet.heat_capacity_cp_j_kg_k") {
    config.planet.heat_capacity_cp_j_kg_k = parse_real(value, line, key);
  } else if (key == "planet.reference_pressure_pa") {
    config.planet.reference_pressure_pa = parse_real(value, line, key);
  } else if (key.starts_with("orbit.") || key == "star.flux_at_semimajor_axis_w_m2") {
    if (!config.orbit.has_value()) config.orbit.emplace();
    if (key == "orbit.period_s")
      config.orbit->period_s = parse_real(value, line, key);
    else if (key == "orbit.eccentricity")
      config.orbit->eccentricity = parse_real(value, line, key);
    else if (key == "orbit.obliquity_rad")
      config.orbit->obliquity_rad = parse_real(value, line, key);
    else if (key == "orbit.longitude_of_periapsis_rad")
      config.orbit->longitude_of_periapsis_rad = parse_real(value, line, key);
    else if (key == "orbit.initial_mean_anomaly_rad")
      config.orbit->initial_mean_anomaly_rad = parse_real(value, line, key);
    else if (key == "orbit.initial_substellar_longitude_rad")
      config.orbit->initial_substellar_longitude_rad = parse_real(value, line, key);
    else
      config.orbit->stellar_flux_at_semimajor_axis_w_m2 = parse_real(value, line, key);
  } else if (key == "run.start_time_s") {
    config.run.start_time_s = parse_real(value, line, key);
  } else if (key == "run.end_time_s") {
    config.run.end_time_s = parse_real(value, line, key);
  } else if (key == "run.time_step_s") {
    config.run.time_step_s = parse_real(value, line, key);
  } else if (key == "run.random_seed") {
    config.run.random_seed = parse_seed(value, line, key);
  } else if (key == "ode.initial_value") {
    config.ode.initial_value = parse_real(value, line, key);
  } else if (key == "ode.decay_rate_s_1") {
    config.ode.decay_rate_s_1 = parse_real(value, line, key);
  } else if (key == "grid.cells_per_panel") {
    config.grid.cells_per_panel = parse_index(value, line, key);
  } else if (key == "transport.test_case") {
    if (value == "solid_body") {
      config.transport.test_case = TransportTestCase::kSolidBody;
    } else if (value == "deformational") {
      config.transport.test_case = TransportTestCase::kDeformational;
    } else if (value == "divergent") {
      config.transport.test_case = TransportTestCase::kDivergent;
    } else {
      throw parse_error(line, "unknown transport.test_case " + std::string(value));
    }
  } else if (key == "transport.initial_condition") {
    if (value == "constant") {
      config.transport.initial_condition = InitialConditionKind::kConstant;
    } else if (value == "gaussian_hill") {
      config.transport.initial_condition = InitialConditionKind::kGaussianHill;
    } else if (value == "cosine_bell") {
      config.transport.initial_condition = InitialConditionKind::kCosineBell;
    } else if (value == "slotted_cylinder") {
      config.transport.initial_condition = InitialConditionKind::kSlottedCylinder;
    } else {
      throw parse_error(line,
                        "unknown transport.initial_condition " + std::string(value));
    }
  } else if (key == "transport.scheme") {
    if (value == "upwind") {
      config.transport.scheme = TransportScheme::kUpwind;
    } else if (value == "linear") {
      config.transport.scheme = TransportScheme::kLinear;
    } else {
      throw parse_error(line, "unknown transport.scheme " + std::string(value));
    }
  } else if (key == "transport.limiter") {
    if (value == "none") {
      config.transport.limiter = LimiterKind::kNone;
    } else if (value == "barth_jespersen") {
      config.transport.limiter = LimiterKind::kBarthJespersen;
    } else {
      throw parse_error(line, "unknown transport.limiter " + std::string(value));
    }
  } else if (key == "transport.cfl") {
    config.transport.cfl = parse_real(value, line, key);
  } else if (key == "transport.rotation_axis_x") {
    config.transport.rotation_axis_x = parse_real(value, line, key);
  } else if (key == "transport.rotation_axis_y") {
    config.transport.rotation_axis_y = parse_real(value, line, key);
  } else if (key == "transport.rotation_axis_z") {
    config.transport.rotation_axis_z = parse_real(value, line, key);
  } else if (key == "transport.angular_speed_rad_s") {
    config.transport.angular_speed_rad_s = parse_real(value, line, key);
  } else if (key == "shallow_water.test_case") {
    if (value == "rest") {
      config.shallow_water.test_case = ShallowWaterTestCase::kRest;
    } else if (value == "linear_wave") {
      config.shallow_water.test_case = ShallowWaterTestCase::kLinearWave;
    } else if (value == "geostrophic_adjustment") {
      config.shallow_water.test_case = ShallowWaterTestCase::kGeostrophicAdjustment;
    } else if (value == "williamson2") {
      config.shallow_water.test_case = ShallowWaterTestCase::kWilliamson2;
    } else if (value == "williamson5") {
      config.shallow_water.test_case = ShallowWaterTestCase::kWilliamson5;
    } else if (value == "williamson6") {
      config.shallow_water.test_case = ShallowWaterTestCase::kWilliamson6;
    } else if (value == "galewsky") {
      config.shallow_water.test_case = ShallowWaterTestCase::kGalewsky;
    } else {
      throw parse_error(line, "unknown shallow_water.test_case " + std::string(value));
    }
  } else if (key == "shallow_water.scheme") {
    if (value == "rusanov") {
      config.shallow_water.scheme = ShallowWaterScheme::kRusanov;
    } else if (value == "compatible") {
      config.shallow_water.scheme = ShallowWaterScheme::kCompatible;
    } else {
      throw parse_error(line, "unknown shallow_water.scheme " + std::string(value));
    }
  } else if (key == "shallow_water.reconstruction") {
    if (value == "piecewise_constant") {
      config.shallow_water.reconstruction = ReconstructionKind::kPiecewiseConstant;
    } else if (value == "linear") {
      config.shallow_water.reconstruction = ReconstructionKind::kLinear;
    } else {
      throw parse_error(line,
                        "unknown shallow_water.reconstruction " + std::string(value));
    }
  } else if (key == "shallow_water.limiter") {
    if (value == "none") {
      config.shallow_water.limiter = LimiterKind::kNone;
    } else if (value == "barth_jespersen") {
      config.shallow_water.limiter = LimiterKind::kBarthJespersen;
    } else {
      throw parse_error(line, "unknown shallow_water.limiter " + std::string(value));
    }
  } else if (key == "shallow_water.cfl") {
    config.shallow_water.cfl = parse_real(value, line, key);
  } else if (key == "shallow_water.mean_depth_m") {
    config.shallow_water.mean_depth_m = parse_real(value, line, key);
  } else if (key == "shallow_water.depth_floor_m") {
    config.shallow_water.depth_floor_m = parse_real(value, line, key);
  } else if (key == "shallow_water.diffusion_kind") {
    if (value == "none") {
      config.shallow_water.diffusion_kind = DiffusionKind::kNone;
    } else if (value == "laplacian") {
      config.shallow_water.diffusion_kind = DiffusionKind::kLaplacian;
    } else if (value == "biharmonic") {
      config.shallow_water.diffusion_kind = DiffusionKind::kBiharmonic;
    } else {
      throw parse_error(line,
                        "unknown shallow_water.diffusion_kind " + std::string(value));
    }
  } else if (key == "shallow_water.diffusion_coefficient") {
    config.shallow_water.diffusion_coefficient = parse_real(value, line, key);
  } else if (key == "shallow_water.flow_axis_x") {
    config.shallow_water.flow_axis_x = parse_real(value, line, key);
  } else if (key == "shallow_water.flow_axis_y") {
    config.shallow_water.flow_axis_y = parse_real(value, line, key);
  } else if (key == "shallow_water.flow_axis_z") {
    config.shallow_water.flow_axis_z = parse_real(value, line, key);
  } else if (key == "shallow_water.maximum_velocity_m_s") {
    config.shallow_water.maximum_velocity_m_s = parse_real(value, line, key);
  } else if (key == "diagnostics.interval_steps") {
    config.diagnostics.interval_steps = parse_seed(value, line, key);
  } else if (key == "vertical.test_case") {
    if (value == "isothermal") {
      config.vertical.test_case = VerticalTestCase::kIsothermal;
    } else if (value == "dry_adiabatic") {
      config.vertical.test_case = VerticalTestCase::kDryAdiabatic;
    } else if (value == "moving_surface_pressure") {
      config.vertical.test_case = VerticalTestCase::kMovingSurfacePressure;
    } else if (value == "manufactured_transport") {
      config.vertical.test_case = VerticalTestCase::kManufacturedTransport;
    } else {
      throw parse_error(line, "unknown vertical.test_case " + std::string(value));
    }
  } else if (key == "vertical.levels") {
    config.vertical.levels = parse_index(value, line, key);
  } else if (key == "vertical.a_half_pa") {
    config.vertical.a_half_pa = parse_real_list(value, line, key);
  } else if (key == "vertical.b_half") {
    config.vertical.b_half = parse_real_list(value, line, key);
  } else if (key == "vertical.surface_pressure_pa") {
    config.vertical.surface_pressure_pa = parse_real(value, line, key);
  } else if (key == "vertical.minimum_surface_pressure_pa") {
    config.vertical.minimum_surface_pressure_pa = parse_real(value, line, key);
  } else if (key == "vertical.maximum_surface_pressure_pa") {
    config.vertical.maximum_surface_pressure_pa = parse_real(value, line, key);
  } else if (key == "vertical.minimum_pressure_thickness_pa") {
    config.vertical.minimum_pressure_thickness_pa = parse_real(value, line, key);
  } else if (key == "vertical.surface_geopotential_m2_s2") {
    config.vertical.surface_geopotential_m2_s2 = parse_real(value, line, key);
  } else if (key == "vertical.initial_temperature_k") {
    config.vertical.initial_temperature_k = parse_real(value, line, key);
  } else if (key == "vertical.initial_potential_temperature_k") {
    config.vertical.initial_potential_temperature_k = parse_real(value, line, key);
  } else if (key == "vertical.temperature_floor_k") {
    config.vertical.temperature_floor_k = parse_real(value, line, key);
  } else if (key == "vertical.transport_scheme") {
    if (value == "donor_cell") {
      config.vertical.transport_scheme = VerticalTransportScheme::kDonorCell;
    } else if (value == "linear") {
      config.vertical.transport_scheme = VerticalTransportScheme::kLinear;
    } else {
      throw parse_error(line,
                        "unknown vertical.transport_scheme " + std::string(value));
    }
  } else if (key == "vertical.limiter") {
    if (value == "none") {
      config.vertical.limiter = VerticalLimiterKind::kNone;
    } else if (value == "minmod") {
      config.vertical.limiter = VerticalLimiterKind::kMinmod;
    } else {
      throw parse_error(line, "unknown vertical.limiter " + std::string(value));
    }
  } else if (key == "vertical.cfl") {
    config.vertical.cfl = parse_real(value, line, key);
  } else if (key == "vertical.forcing_amplitude") {
    config.vertical.forcing_amplitude = parse_real(value, line, key);
  } else if (key == "dry_hydrostatic.test_case") {
    if (value == "isothermal_rest")
      config.dry_hydrostatic.test_case = DryHydrostaticTestCase::kIsothermalRest;
    else if (value == "solid_body_transport")
      config.dry_hydrostatic.test_case = DryHydrostaticTestCase::kSolidBodyTransport;
    else if (value == "dcmip_deformational")
      config.dry_hydrostatic.test_case = DryHydrostaticTestCase::kDcmipDeformational;
    else if (value == "dcmip_hadley")
      config.dry_hydrostatic.test_case = DryHydrostaticTestCase::kDcmipHadley;
    else if (value == "linear_wave")
      config.dry_hydrostatic.test_case = DryHydrostaticTestCase::kLinearWave;
    else if (value == "dcmip_2_0_0_rest")
      config.dry_hydrostatic.test_case = DryHydrostaticTestCase::kDcmip200Rest;
    else if (value == "linear_mountain_wave")
      config.dry_hydrostatic.test_case = DryHydrostaticTestCase::kLinearMountainWave;
    else if (value == "jw06_steady")
      config.dry_hydrostatic.test_case = DryHydrostaticTestCase::kJw06Steady;
    else if (value == "jw06_baroclinic")
      config.dry_hydrostatic.test_case = DryHydrostaticTestCase::kJw06Baroclinic;
    else if (value == "umjs14_steady")
      config.dry_hydrostatic.test_case = DryHydrostaticTestCase::kUmjs14Steady;
    else if (value == "umjs14_baroclinic")
      config.dry_hydrostatic.test_case = DryHydrostaticTestCase::kUmjs14Baroclinic;
    else if (value == "held_suarez")
      config.dry_hydrostatic.test_case = DryHydrostaticTestCase::kHeldSuarez;
    else
      throw parse_error(line,
                        "unknown dry_hydrostatic.test_case " + std::string(value));
  } else if (key == "dry_hydrostatic.reconstruction") {
    if (value == "linear")
      config.dry_hydrostatic.reconstruction = ReconstructionKind::kLinear;
    else if (value == "piecewise_constant")
      config.dry_hydrostatic.reconstruction = ReconstructionKind::kPiecewiseConstant;
    else
      throw parse_error(line,
                        "unknown dry_hydrostatic.reconstruction " + std::string(value));
  } else if (key == "dry_hydrostatic.limiter") {
    if (value == "barth_jespersen")
      config.dry_hydrostatic.limiter = LimiterKind::kBarthJespersen;
    else if (value == "none")
      config.dry_hydrostatic.limiter = LimiterKind::kNone;
    else
      throw parse_error(line, "unknown dry_hydrostatic.limiter " + std::string(value));
  } else if (key == "dry_hydrostatic.cfl") {
    config.dry_hydrostatic.cfl = parse_real(value, line, key);
  } else if (key == "dry_hydrostatic.diffusion_kind") {
    if (value == "none")
      config.dry_hydrostatic.diffusion_kind = DiffusionKind::kNone;
    else if (value == "laplacian")
      config.dry_hydrostatic.diffusion_kind = DiffusionKind::kLaplacian;
    else if (value == "biharmonic")
      config.dry_hydrostatic.diffusion_kind = DiffusionKind::kBiharmonic;
    else
      throw parse_error(line,
                        "unknown dry_hydrostatic.diffusion_kind " + std::string(value));
  } else if (key == "dry_hydrostatic.diffusion_coefficient") {
    config.dry_hydrostatic.diffusion_coefficient = parse_real(value, line, key);
  } else if (key == "dry_hydrostatic.time_integrator") {
    if (value != "semi_implicit")
      throw parse_error(
          line, "unknown dry_hydrostatic.time_integrator " + std::string(value));
    config.dry_hydrostatic.time_integrator =
        DryHydrostaticTimeIntegrator::kSemiImplicit;
  } else if (key == "dry_hydrostatic.advective_cfl") {
    config.dry_hydrostatic.advective_cfl = parse_real(value, line, key);
  } else if (key.starts_with("semi_implicit.")) {
    if (!config.semi_implicit.has_value()) config.semi_implicit.emplace();
    auto& semi_implicit = *config.semi_implicit;
    if (key == "semi_implicit.reference_surface_pressure_pa")
      semi_implicit.reference_surface_pressure_pa = parse_real(value, line, key);
    else if (key == "semi_implicit.reference_temperature_k")
      semi_implicit.reference_temperature_k = parse_real(value, line, key);
    else if (key == "semi_implicit.reference_update") {
      if (value == "fixed")
        semi_implicit.reference_update = SemiImplicitReferenceUpdate::kFixed;
      else if (value == "per_step")
        semi_implicit.reference_update = SemiImplicitReferenceUpdate::kPerStep;
      else
        throw parse_error(
            line, "unknown semi_implicit.reference_update " + std::string(value));
    } else if (key == "semi_implicit.implicit_weight")
      semi_implicit.implicit_weight = parse_real(value, line, key);
    else if (key == "semi_implicit.wave_cfl_threshold")
      semi_implicit.wave_cfl_threshold = parse_real(value, line, key);
    else if (key == "semi_implicit.maximum_implicit_modes")
      semi_implicit.maximum_implicit_modes = parse_index(value, line, key);
    else if (key == "semi_implicit.nonlinear_iterations")
      semi_implicit.nonlinear_iterations = parse_index(value, line, key);
    else if (key == "semi_implicit.linear_relative_tolerance")
      semi_implicit.linear_relative_tolerance = parse_real(value, line, key);
    else if (key == "semi_implicit.linear_absolute_tolerance")
      semi_implicit.linear_absolute_tolerance = parse_real(value, line, key);
    else if (key == "semi_implicit.linear_maximum_iterations")
      semi_implicit.linear_maximum_iterations = parse_index(value, line, key);
    else if (key == "semi_implicit.gmres_restart")
      semi_implicit.gmres_restart = parse_index(value, line, key);
    else if (key == "semi_implicit.minimum_time_step_s")
      semi_implicit.minimum_time_step_s = parse_real(value, line, key);
    else
      throw parse_error(line, "unknown key " + std::string(key));
  } else if (key == "orography.kind") {
    if (value == "dcmip_2_0_0")
      config.orography.kind = OrographyKind::kDcmip200;
    else if (value == "williamson5")
      config.orography.kind = OrographyKind::kWilliamson5;
    else if (value == "linear_bell")
      config.orography.kind = OrographyKind::kLinearBell;
    else if (value == "jw06")
      config.orography.kind = OrographyKind::kJw06;
    else if (value == "latlon_csv")
      config.orography.kind = OrographyKind::kLatLonCsv;
    else
      throw parse_error(line, "unknown orography.kind " + std::string(value));
  } else if (key == "physics.kind") {
    if (value == "none")
      config.physics.kind = PhysicsKind::kNone;
    else if (value == "held_suarez")
      config.physics.kind = PhysicsKind::kHeldSuarez;
    else if (value == "planetary_newtonian")
      config.physics.kind = PhysicsKind::kPlanetaryNewtonian;
    else if (value == "surface_energy_balance")
      config.physics.kind = PhysicsKind::kSurfaceEnergyBalance;
    else if (value == "gray_radiation")
      config.physics.kind = PhysicsKind::kGrayRadiation;
    else
      throw parse_error(line, "unknown physics.kind " + std::string(value));
  } else if (key == "forcing.geometry") {
    if (value == "axisymmetric")
      config.physics.geometry = ForcingGeometry::kAxisymmetric;
    else if (value == "substellar")
      config.physics.geometry = ForcingGeometry::kSubstellar;
    else
      throw parse_error(line, "unknown forcing.geometry " + std::string(value));
  } else if (key == "surface.geography") {
    if (!config.surface.has_value()) config.surface.emplace();
    auto& surface = *config.surface;
    if (value == "uniform")
      surface.geography = SurfaceGeography::kUniform;
    else if (value == "earth")
      surface.geography = SurfaceGeography::kEarth;
    else
      throw parse_error(line, "unknown surface.geography " + std::string(value));
  } else if (key == "surface.uniform_land_fraction") {
    if (!config.surface.has_value()) config.surface.emplace();
    config.surface->uniform_land_fraction = parse_real(value, line, key);
  } else if (key == "surface.input_file") {
    if (!config.surface.has_value()) config.surface.emplace();
    config.surface->input_file = value;
  } else if (key == "surface.input_fingerprint_fnv1a64") {
    if (!config.surface.has_value()) config.surface.emplace();
    config.surface->input_fingerprint_fnv1a64 = value;
  } else if (key == "surface.quadrature_order") {
    if (!config.surface.has_value()) config.surface.emplace();
    config.surface->quadrature_order = parse_index(value, line, key);
  } else if (key == "surface.smoothing_passes") {
    if (!config.surface.has_value()) config.surface.emplace();
    config.surface->smoothing_passes = parse_index(value, line, key);
  } else if (key.starts_with("surface.")) {
    if (!config.surface.has_value()) config.surface.emplace();
    const Real parsed = parse_real(value, line, key);
    if (key == "surface.land_heat_capacity_j_m2_k")
      config.surface->land_heat_capacity_j_m2_k = parsed;
    else if (key == "surface.ocean_heat_capacity_j_m2_k")
      config.surface->ocean_heat_capacity_j_m2_k = parsed;
    else if (key == "surface.initial_temperature_k")
      config.surface->initial_temperature_k = parsed;
    else if (key == "surface.albedo")
      config.surface->albedo = parsed;
    else if (key == "surface.emissivity")
      config.surface->emissivity = parsed;
    else if (key == "surface.air_exchange_coefficient_w_m2_k")
      config.surface->air_exchange_coefficient_w_m2_k = parsed;
    else if (key == "surface.internal_heat_flux_w_m2")
      config.surface->internal_heat_flux_w_m2 = parsed;
    else if (key == "surface.cfl")
      config.surface->cfl = parsed;
    else
      throw parse_error(line, "unknown key " + std::string(key));
  } else if (key.starts_with("radiation.")) {
    if (!config.radiation.has_value()) config.radiation.emplace();
    const Real parsed = parse_real(value, line, key);
    if (key == "radiation.shortwave_absorption_m2_kg")
      config.radiation->shortwave_absorption_m2_kg = parsed;
    else if (key == "radiation.longwave_absorption_ref_m2_kg")
      config.radiation->longwave_absorption_ref_m2_kg = parsed;
    else if (key == "radiation.reference_pressure_pa")
      config.radiation->reference_pressure_pa = parsed;
    else if (key == "radiation.longwave_pressure_exponent")
      config.radiation->longwave_pressure_exponent = parsed;
    else if (key == "radiation.longwave_diffusivity_factor")
      config.radiation->longwave_diffusivity_factor = parsed;
    else if (key == "radiation.shortwave_diffuse_factor")
      config.radiation->shortwave_diffuse_factor = parsed;
    else if (key == "radiation.cfl")
      config.radiation->cfl = parsed;
    else
      throw parse_error(line, "unknown key " + std::string(key));
  } else if (key == "orography.input_file") {
    config.orography.input_file = value;
  } else if (key == "orography.input_fingerprint_fnv1a64") {
    config.orography.input_fingerprint_fnv1a64 = value;
  } else if (key == "orography.smoothing_passes") {
    config.orography.smoothing_passes = parse_index(value, line, key);
  } else if (key == "output.directory") {
    config.output_directory = value;
  } else {
    throw parse_error(line, "unknown key " + std::string(key));
  }
}

template <std::size_t Size>
void require_keys(const std::set<std::string, std::less<>>& seen_keys,
                  const std::array<std::string_view, Size>& keys) {
  for (const auto key : keys) {
    if (!seen_keys.contains(key)) {
      throw std::runtime_error("configuration is missing required key " +
                               std::string(key));
    }
  }
}

}  // namespace

void ExperimentConfig::validate() const {
  planet.validate();
  if (orbit.has_value()) orbit->validate();
  if (kind != ExperimentKind::kDryHydrostatic &&
      (dry_hydrostatic.time_integrator == DryHydrostaticTimeIntegrator::kSemiImplicit ||
       semi_implicit.has_value()))
    throw std::invalid_argument(
        "semi-implicit integration is valid only for dry_hydrostatic");
  require_non_negative(run.start_time_s, "run.start_time_s");
  require_finite(run.end_time_s, "run.end_time_s");
  require_positive(run.time_step_s, "run.time_step_s");
  if (kind == ExperimentKind::kOde) {
    require_finite(ode.initial_value, "ode.initial_value");
    require_non_negative(ode.decay_rate_s_1, "ode.decay_rate_s_1");
  } else if (kind == ExperimentKind::kSphereTransport) {
    if (grid.cells_per_panel <= 0) {
      throw std::invalid_argument("grid.cells_per_panel must be positive");
    }
    require_finite(transport.cfl, "transport.cfl");
    if (!(transport.cfl > 0.0 && transport.cfl <= 1.0)) {
      throw std::invalid_argument("transport.cfl must be in (0, 1]");
    }
    require_finite(transport.rotation_axis_x, "transport.rotation_axis_x");
    require_finite(transport.rotation_axis_y, "transport.rotation_axis_y");
    require_finite(transport.rotation_axis_z, "transport.rotation_axis_z");
    require_finite(transport.angular_speed_rad_s, "transport.angular_speed_rad_s");
    const Real axis_norm_squared =
        transport.rotation_axis_x * transport.rotation_axis_x +
        transport.rotation_axis_y * transport.rotation_axis_y +
        transport.rotation_axis_z * transport.rotation_axis_z;
    if (!(axis_norm_squared > 0.0) || !std::isfinite(axis_norm_squared)) {
      throw std::invalid_argument("transport rotation axis must be nonzero");
    }
  } else if (kind == ExperimentKind::kShallowWater) {
    if (grid.cells_per_panel <= 0) {
      throw std::invalid_argument("grid.cells_per_panel must be positive");
    }
    require_finite(shallow_water.cfl, "shallow_water.cfl");
    if (!(shallow_water.cfl > 0.0 && shallow_water.cfl <= 1.0)) {
      throw std::invalid_argument("shallow_water.cfl must be in (0, 1]");
    }
    require_positive(shallow_water.mean_depth_m, "shallow_water.mean_depth_m");
    require_non_negative(shallow_water.depth_floor_m, "shallow_water.depth_floor_m");
    if (!(shallow_water.mean_depth_m > shallow_water.depth_floor_m)) {
      throw std::invalid_argument(
          "shallow_water.mean_depth_m must exceed depth_floor_m");
    }
    require_non_negative(shallow_water.diffusion_coefficient,
                         "shallow_water.diffusion_coefficient");
    if ((shallow_water.diffusion_kind == DiffusionKind::kNone) !=
        (shallow_water.diffusion_coefficient == 0.0)) {
      throw std::invalid_argument(
          "shallow_water diffusion coefficient must be zero exactly when kind is none");
    }
    require_finite(shallow_water.flow_axis_x, "shallow_water.flow_axis_x");
    require_finite(shallow_water.flow_axis_y, "shallow_water.flow_axis_y");
    require_finite(shallow_water.flow_axis_z, "shallow_water.flow_axis_z");
    const Real axis_norm_squared =
        shallow_water.flow_axis_x * shallow_water.flow_axis_x +
        shallow_water.flow_axis_y * shallow_water.flow_axis_y +
        shallow_water.flow_axis_z * shallow_water.flow_axis_z;
    if (!(axis_norm_squared > 0.0) || !std::isfinite(axis_norm_squared)) {
      throw std::invalid_argument("shallow-water flow axis must be nonzero");
    }
    require_non_negative(shallow_water.maximum_velocity_m_s,
                         "shallow_water.maximum_velocity_m_s");
    if (diagnostics.interval_steps == 0) {
      throw std::invalid_argument("diagnostics.interval_steps must be positive");
    }
    if ((shallow_water.test_case == ShallowWaterTestCase::kWilliamson5) !=
        (orography.kind == OrographyKind::kWilliamson5))
      throw std::invalid_argument(
          "williamson5 test case requires williamson5 orography and vice versa");
    if (orography.kind != OrographyKind::kFlat &&
        shallow_water.scheme == ShallowWaterScheme::kCompatible)
      throw std::invalid_argument(
          "non-flat shallow-water runs require the Rusanov scheme");
  } else {
    if (vertical.levels <= 0) {
      throw std::invalid_argument("vertical.levels must be positive");
    }
    const auto expected_size = static_cast<std::size_t>(vertical.levels + 1);
    if (vertical.a_half_pa.size() != expected_size ||
        vertical.b_half.size() != expected_size) {
      throw std::invalid_argument("vertical A/B lists must contain levels + 1 values");
    }
    for (const auto value : vertical.a_half_pa) {
      require_finite(value, "vertical.a_half_pa");
    }
    for (const auto value : vertical.b_half) {
      require_finite(value, "vertical.b_half");
      if (value < 0.0 || value > 1.0) {
        throw std::invalid_argument("vertical.b_half must be in [0, 1]");
      }
    }
    require_positive(vertical.surface_pressure_pa, "vertical.surface_pressure_pa");
    require_positive(vertical.minimum_surface_pressure_pa,
                     "vertical.minimum_surface_pressure_pa");
    require_positive(vertical.maximum_surface_pressure_pa,
                     "vertical.maximum_surface_pressure_pa");
    if (vertical.minimum_surface_pressure_pa > vertical.surface_pressure_pa ||
        vertical.surface_pressure_pa > vertical.maximum_surface_pressure_pa) {
      throw std::invalid_argument("vertical.surface_pressure_pa is outside its range");
    }
    require_positive(vertical.minimum_pressure_thickness_pa,
                     "vertical.minimum_pressure_thickness_pa");
    require_finite(vertical.surface_geopotential_m2_s2,
                   "vertical.surface_geopotential_m2_s2");
    require_positive(vertical.initial_temperature_k, "vertical.initial_temperature_k");
    require_positive(vertical.initial_potential_temperature_k,
                     "vertical.initial_potential_temperature_k");
    require_positive(vertical.temperature_floor_k, "vertical.temperature_floor_k");
    require_finite(vertical.cfl, "vertical.cfl");
    if (!(vertical.cfl > 0.0 && vertical.cfl <= 1.0)) {
      throw std::invalid_argument("vertical.cfl must be in (0, 1]");
    }
    require_finite(vertical.forcing_amplitude, "vertical.forcing_amplitude");
    if (vertical.a_half_pa.front() <= 0.0 ||
        vertical.minimum_surface_pressure_pa <= vertical.a_half_pa.front()) {
      throw std::invalid_argument(
          "vertical model top must be positive and below ps range");
    }
    if (vertical.b_half.front() != 0.0 || vertical.a_half_pa.back() != 0.0 ||
        vertical.b_half.back() != 1.0) {
      throw std::invalid_argument(
          "vertical hybrid endpoints violate atmospheric convention");
    }
    for (std::size_t index = 1; index < expected_size; ++index) {
      if (vertical.b_half[index] < vertical.b_half[index - 1]) {
        throw std::invalid_argument("vertical.b_half must be nondecreasing");
      }
      for (const auto surface_pressure : {vertical.minimum_surface_pressure_pa,
                                          vertical.maximum_surface_pressure_pa}) {
        const auto thickness =
            (vertical.a_half_pa[index] + vertical.b_half[index] * surface_pressure) -
            (vertical.a_half_pa[index - 1] +
             vertical.b_half[index - 1] * surface_pressure);
        if (!(thickness >= vertical.minimum_pressure_thickness_pa)) {
          throw std::invalid_argument("vertical layer is nonmonotone or too thin");
        }
      }
    }
    if (diagnostics.interval_steps == 0) {
      throw std::invalid_argument("diagnostics.interval_steps must be positive");
    }
  }

  if (kind == ExperimentKind::kDryHydrostatic) {
    if (grid.cells_per_panel <= 0)
      throw std::invalid_argument("grid.cells_per_panel must be positive");
    if (vertical.surface_geopotential_m2_s2 != 0.0)
      throw std::invalid_argument("dry_hydrostatic requires flat surface geopotential");
    require_finite(dry_hydrostatic.cfl, "dry_hydrostatic.cfl");
    if (!(dry_hydrostatic.cfl > 0.0 && dry_hydrostatic.cfl <= 1.0))
      throw std::invalid_argument("dry_hydrostatic.cfl must be in (0, 1]");
    const bool uses_semi_implicit =
        dry_hydrostatic.time_integrator == DryHydrostaticTimeIntegrator::kSemiImplicit;
    if (uses_semi_implicit != semi_implicit.has_value())
      throw std::invalid_argument(
          "semi-implicit integrator and parameters must be configured together");
    if (uses_semi_implicit) {
      const auto& parameters = *semi_implicit;
      require_finite(dry_hydrostatic.advective_cfl, "dry_hydrostatic.advective_cfl");
      if (!(dry_hydrostatic.advective_cfl > 0.0 &&
            dry_hydrostatic.advective_cfl <= 1.0))
        throw std::invalid_argument("dry_hydrostatic.advective_cfl must be in (0, 1]");
      require_positive(parameters.reference_surface_pressure_pa,
                       "semi_implicit.reference_surface_pressure_pa");
      if (parameters.reference_surface_pressure_pa <
              vertical.minimum_surface_pressure_pa ||
          parameters.reference_surface_pressure_pa >
              vertical.maximum_surface_pressure_pa)
        throw std::invalid_argument(
            "semi-implicit reference surface pressure is outside configured bounds");
      require_positive(parameters.reference_temperature_k,
                       "semi_implicit.reference_temperature_k");
      if (parameters.reference_temperature_k < vertical.temperature_floor_k)
        throw std::invalid_argument(
            "semi-implicit reference temperature is below the temperature floor");
      require_finite(parameters.implicit_weight, "semi_implicit.implicit_weight");
      if (parameters.implicit_weight < 0.5 || parameters.implicit_weight > 0.55)
        throw std::invalid_argument(
            "semi_implicit.implicit_weight must be in [0.5, 0.55]");
      require_finite(parameters.wave_cfl_threshold, "semi_implicit.wave_cfl_threshold");
      if (!(parameters.wave_cfl_threshold > 0.0 &&
            parameters.wave_cfl_threshold <= 1.0))
        throw std::invalid_argument(
            "semi_implicit.wave_cfl_threshold must be in (0, 1]");
      if (parameters.maximum_implicit_modes <= 0 ||
          parameters.maximum_implicit_modes > vertical.levels)
        throw std::invalid_argument(
            "semi_implicit.maximum_implicit_modes must be in [1, vertical.levels]");
      // Benard (2003): a non-extrapolating ICI scheme with a single iteration is only
      // first-order accurate in time, so two is the smallest second-order setting.
      if (parameters.nonlinear_iterations < 2)
        throw std::invalid_argument(
            "semi_implicit.nonlinear_iterations must be at least 2");
      require_positive(parameters.linear_relative_tolerance,
                       "semi_implicit.linear_relative_tolerance");
      require_positive(parameters.linear_absolute_tolerance,
                       "semi_implicit.linear_absolute_tolerance");
      if (parameters.linear_maximum_iterations <= 0)
        throw std::invalid_argument(
            "semi_implicit.linear_maximum_iterations must be positive");
      if (parameters.gmres_restart <= 0 ||
          parameters.gmres_restart > parameters.linear_maximum_iterations)
        throw std::invalid_argument(
            "semi_implicit.gmres_restart must be in [1, linear maximum]");
      require_positive(parameters.minimum_time_step_s,
                       "semi_implicit.minimum_time_step_s");
      if (parameters.minimum_time_step_s > run.time_step_s)
        throw std::invalid_argument(
            "semi_implicit.minimum_time_step_s must not exceed run.time_step_s");
    }
    require_non_negative(dry_hydrostatic.diffusion_coefficient,
                         "dry_hydrostatic.diffusion_coefficient");
    if ((dry_hydrostatic.diffusion_kind == DiffusionKind::kNone) !=
        (dry_hydrostatic.diffusion_coefficient == 0.0))
      throw std::invalid_argument(
          "dry_hydrostatic diffusion coefficient must be zero exactly when kind is "
          "none");
    if ((dry_hydrostatic.test_case == DryHydrostaticTestCase::kDcmip200Rest) !=
        (orography.kind == OrographyKind::kDcmip200))
      throw std::invalid_argument(
          "dcmip_2_0_0_rest requires dcmip_2_0_0 orography and vice versa");
    if ((dry_hydrostatic.test_case == DryHydrostaticTestCase::kLinearMountainWave) !=
        (orography.kind == OrographyKind::kLinearBell))
      throw std::invalid_argument(
          "linear_mountain_wave requires linear_bell orography and vice versa");
    const bool jw06_case =
        dry_hydrostatic.test_case == DryHydrostaticTestCase::kJw06Steady ||
        dry_hydrostatic.test_case == DryHydrostaticTestCase::kJw06Baroclinic;
    if (jw06_case != (orography.kind == OrographyKind::kJw06))
      throw std::invalid_argument(
          "JW06 test cases require jw06 orography and vice versa");
    if (jw06_case &&
        (planet.radius_m != kJw06Planet.radius_m ||
         planet.rotation_rate_rad_s != kJw06Planet.rotation_rate_rad_s ||
         planet.gravity_m_s2 != kJw06Planet.gravity_m_s2 ||
         planet.gas_constant_j_kg_k != kJw06Planet.gas_constant_j_kg_k ||
         planet.heat_capacity_cp_j_kg_k != kJw06Planet.heat_capacity_cp_j_kg_k ||
         planet.reference_pressure_pa != kJw06Planet.reference_pressure_pa))
      throw std::invalid_argument("JW06 requires registered Earth constants");
    const bool held_suarez_case =
        dry_hydrostatic.test_case == DryHydrostaticTestCase::kHeldSuarez;
    const bool benchmark_physics = physics.kind == PhysicsKind::kHeldSuarez ||
                                   physics.kind == PhysicsKind::kPlanetaryNewtonian;
    if (held_suarez_case != benchmark_physics)
      throw std::invalid_argument(
          "held_suarez test case requires analytic benchmark physics and vice versa");
    if (physics.kind == PhysicsKind::kHeldSuarez) {
      if (orography.kind != OrographyKind::kFlat)
        throw std::invalid_argument("held_suarez requires flat orography");
      constexpr PlanetParameters earth{6371220.0, 7.29212e-5, 9.80616,
                                       287.0,     1004.0,     100000.0};
      if (planet.radius_m != earth.radius_m ||
          planet.rotation_rate_rad_s != earth.rotation_rate_rad_s ||
          planet.gravity_m_s2 != earth.gravity_m_s2 ||
          planet.gas_constant_j_kg_k != earth.gas_constant_j_kg_k ||
          planet.heat_capacity_cp_j_kg_k != earth.heat_capacity_cp_j_kg_k ||
          planet.reference_pressure_pa != earth.reference_pressure_pa)
        throw std::invalid_argument("held_suarez requires registered Earth constants");
    }
    if (physics.kind == PhysicsKind::kPlanetaryNewtonian) {
      if (orography.kind != OrographyKind::kFlat)
        throw std::invalid_argument("planetary_newtonian requires flat orography");
      if (!surface.has_value() || surface->geography != SurfaceGeography::kUniform)
        throw std::invalid_argument(
            "planetary_newtonian requires uniform surface geography");
      if (physics.geometry == ForcingGeometry::kSubstellar && !orbit.has_value())
        throw std::invalid_argument("substellar forcing requires orbit parameters");
    }
    if (physics.kind == PhysicsKind::kSurfaceEnergyBalance ||
        physics.kind == PhysicsKind::kGrayRadiation) {
      if (!surface.has_value() || !orbit.has_value())
        throw std::invalid_argument(
            "surface radiation physics requires surface and orbit parameters");
      require_positive(surface->land_heat_capacity_j_m2_k,
                       "surface.land_heat_capacity_j_m2_k");
      require_positive(surface->ocean_heat_capacity_j_m2_k,
                       "surface.ocean_heat_capacity_j_m2_k");
      require_positive(surface->initial_temperature_k, "surface.initial_temperature_k");
      require_finite(surface->albedo, "surface.albedo");
      require_finite(surface->emissivity, "surface.emissivity");
      if (surface->albedo < 0.0 || surface->albedo > 1.0 || surface->emissivity < 0.0 ||
          surface->emissivity > 1.0)
        throw std::invalid_argument("surface albedo and emissivity must be in [0, 1]");
      require_non_negative(surface->air_exchange_coefficient_w_m2_k,
                           "surface.air_exchange_coefficient_w_m2_k");
      require_finite(surface->internal_heat_flux_w_m2,
                     "surface.internal_heat_flux_w_m2");
      if (!(surface->cfl > 0.0) || !(surface->cfl <= 1.0))
        throw std::invalid_argument("surface.cfl must be in (0, 1]");
    }
    if (physics.kind == PhysicsKind::kGrayRadiation) {
      if (!radiation.has_value())
        throw std::invalid_argument("gray_radiation requires radiation parameters");
      require_non_negative(radiation->shortwave_absorption_m2_kg,
                           "radiation.shortwave_absorption_m2_kg");
      require_non_negative(radiation->longwave_absorption_ref_m2_kg,
                           "radiation.longwave_absorption_ref_m2_kg");
      require_positive(radiation->reference_pressure_pa,
                       "radiation.reference_pressure_pa");
      require_finite(radiation->longwave_pressure_exponent,
                     "radiation.longwave_pressure_exponent");
      if (radiation->longwave_pressure_exponent < 1.0)
        throw std::invalid_argument(
            "radiation.longwave_pressure_exponent must be at least one");
      require_positive(radiation->longwave_diffusivity_factor,
                       "radiation.longwave_diffusivity_factor");
      require_positive(radiation->shortwave_diffuse_factor,
                       "radiation.shortwave_diffuse_factor");
      require_finite(radiation->cfl, "radiation.cfl");
      if (!(radiation->cfl > 0.0 && radiation->cfl <= 1.0))
        throw std::invalid_argument("radiation.cfl must be in (0, 1]");
    } else if (radiation.has_value()) {
      throw std::invalid_argument(
          "radiation parameters are valid only for gray_radiation");
    }
    if (surface.has_value()) {
      require_finite(surface->uniform_land_fraction, "surface.uniform_land_fraction");
      if (surface->uniform_land_fraction < 0.0 || surface->uniform_land_fraction > 1.0)
        throw std::invalid_argument("surface.uniform_land_fraction must be in [0, 1]");
      if (surface->geography == SurfaceGeography::kEarth &&
          surface->uniform_land_fraction != 0.0)
        throw std::invalid_argument(
            "earth geography does not accept a uniform land fraction");
      if (surface->geography == SurfaceGeography::kUniform &&
          (!surface->input_file.empty() ||
           !surface->input_fingerprint_fnv1a64.empty() ||
           surface->quadrature_order != 1 || surface->smoothing_passes != 0))
        throw std::invalid_argument(
            "uniform surface geography does not accept data options");
      if (surface->geography == SurfaceGeography::kEarth) {
        if (surface->input_file.empty() ||
            surface->input_file.find_first_of("\r\n") != std::string::npos ||
            std::filesystem::path(surface->input_file).is_absolute())
          throw std::invalid_argument(
              "earth surface geography requires a relative input file");
        if (surface->input_fingerprint_fnv1a64.size() != 16 ||
            !std::ranges::all_of(surface->input_fingerprint_fnv1a64,
                                 [](const char value) {
                                   return (value >= '0' && value <= '9') ||
                                          (value >= 'a' && value <= 'f');
                                 }))
          throw std::invalid_argument(
              "earth surface fingerprint must be 16 lowercase hexadecimal digits");
        if (surface->quadrature_order <= 0 || surface->quadrature_order > 8)
          throw std::invalid_argument("surface.quadrature_order must be in [1, 8]");
        if (surface->smoothing_passes < 0)
          throw std::invalid_argument("surface.smoothing_passes must be nonnegative");
      }
    }
  } else if (physics.kind != PhysicsKind::kNone) {
    throw std::invalid_argument("physics is supported only for dry_hydrostatic");
  } else if (surface.has_value() || radiation.has_value()) {
    throw std::invalid_argument(
        "surface and radiation are supported only for dry_hydrostatic");
  }

  if (orography.kind == OrographyKind::kFlat) {
    if (!orography.input_file.empty() || !orography.input_fingerprint_fnv1a64.empty() ||
        orography.smoothing_passes != 0)
      throw std::invalid_argument("flat orography does not accept input options");
  } else if (orography.kind == OrographyKind::kLatLonCsv) {
    if (orography.input_file.empty() ||
        orography.input_file.find_first_of("\r\n") != std::string::npos ||
        std::filesystem::path(orography.input_file).is_absolute())
      throw std::invalid_argument("latlon_csv requires a single-line input file");
    if (orography.input_fingerprint_fnv1a64.size() != 16 ||
        !std::ranges::all_of(orography.input_fingerprint_fnv1a64, [](char value) {
          return (value >= '0' && value <= '9') || (value >= 'a' && value <= 'f');
        }))
      throw std::invalid_argument(
          "latlon_csv fingerprint must be 16 lowercase hexadecimal digits");
    if (orography.smoothing_passes < 0)
      throw std::invalid_argument("orography smoothing passes must be nonnegative");
  } else if (!orography.input_file.empty() ||
             !orography.input_fingerprint_fnv1a64.empty() ||
             orography.smoothing_passes != 0) {
    throw std::invalid_argument("analytic orography does not accept input options");
  }

  if (run.end_time_s <= run.start_time_s) {
    throw std::invalid_argument("run.end_time_s must be greater than run.start_time_s");
  }
  if (output_directory.empty()) {
    throw std::invalid_argument("output.directory must not be empty");
  }
  if (output_directory.find_first_of("\r\n") != std::string::npos) {
    throw std::invalid_argument("output.directory must be a single line");
  }
}

ExperimentConfig parse_experiment_config(std::istream& input) {
  ExperimentConfig config{};
  std::set<std::string, std::less<>> seen_keys;
  std::string line_text;
  std::size_t line_number = 0;

  while (std::getline(input, line_text)) {
    ++line_number;
    const auto line = trim(line_text);
    if (line.empty() || line.front() == '#') {
      continue;
    }

    const auto separator = line.find('=');
    if (separator == std::string_view::npos ||
        line.find('=', separator + 1) != std::string_view::npos) {
      throw parse_error(line_number, "expected exactly one '=' separator");
    }

    const auto key = trim(line.substr(0, separator));
    const auto value = trim(line.substr(separator + 1));
    if (key.empty() || value.empty()) {
      throw parse_error(line_number, "key and value must not be empty");
    }
    if (!seen_keys.emplace(key).second) {
      throw parse_error(line_number, "duplicate key " + std::string(key));
    }
    assign_value(config, key, value, line_number);
  }

  if (input.bad()) {
    throw std::runtime_error("failed while reading configuration stream");
  }

  require_keys(seen_keys, kCommonRequiredKeys);
  const bool has_orbit_key = std::ranges::any_of(
      kOrbitRequiredKeys, [&](const auto key) { return seen_keys.contains(key); });
  if (has_orbit_key) require_keys(seen_keys, kOrbitRequiredKeys);
  if (config.kind != ExperimentKind::kDryHydrostatic &&
      (has_orbit_key || seen_keys.contains("forcing.geometry")))
    throw std::runtime_error(
        "orbit and forcing keys are valid only for dry_hydrostatic");
  if (config.kind == ExperimentKind::kOde) {
    require_keys(seen_keys, kOdeRequiredKeys);
    for (const auto& key : seen_keys) {
      if (key.starts_with("grid.") || key.starts_with("transport.") ||
          key.starts_with("shallow_water.") || key.starts_with("diagnostics.") ||
          key.starts_with("orography.") || key.starts_with("physics.") ||
          key.starts_with("surface.") || key.starts_with("radiation.")) {
        throw std::runtime_error("key " + key + " is not valid for ode experiment");
      }
    }
  } else if (config.kind == ExperimentKind::kSphereTransport) {
    require_keys(seen_keys, kTransportRequiredKeys);
    for (const auto& key : seen_keys) {
      if (key.starts_with("ode.") || key.starts_with("shallow_water.") ||
          key.starts_with("diagnostics.") || key.starts_with("orography.") ||
          key.starts_with("physics.") || key.starts_with("surface.") ||
          key.starts_with("radiation.")) {
        throw std::runtime_error("key " + key +
                                 " is not valid for sphere_transport experiment");
      }
    }
  } else if (config.kind == ExperimentKind::kShallowWater) {
    require_keys(seen_keys, kShallowWaterRequiredKeys);
    for (const auto& key : seen_keys) {
      if (key.starts_with("ode.") || key.starts_with("transport.") ||
          key.starts_with("physics.") || key.starts_with("surface.") ||
          key.starts_with("radiation.")) {
        throw std::runtime_error("key " + key +
                                 " is not valid for shallow_water experiment");
      }
    }
  } else if (config.kind == ExperimentKind::kVerticalColumn) {
    require_keys(seen_keys, kVerticalRequiredKeys);
    for (const auto& key : seen_keys) {
      if (key.starts_with("ode.") || key.starts_with("grid.") ||
          key.starts_with("transport.") || key.starts_with("shallow_water.") ||
          key.starts_with("orography.") || key.starts_with("physics.") ||
          key.starts_with("surface.") || key.starts_with("radiation.")) {
        throw std::runtime_error("key " + key +
                                 " is not valid for vertical_column experiment");
      }
    }
  } else {
    require_keys(seen_keys, kDryHydrostaticRequiredKeys);
    for (const auto& key : seen_keys) {
      if (key.starts_with("ode.") || key.starts_with("transport.") ||
          key.starts_with("shallow_water.") || key == "vertical.test_case" ||
          key == "vertical.forcing_amplitude") {
        throw std::runtime_error("key " + key +
                                 " is not valid for dry_hydrostatic experiment");
      }
    }
  }

  const bool has_time_integrator =
      seen_keys.contains("dry_hydrostatic.time_integrator");
  const bool has_advective_cfl = seen_keys.contains("dry_hydrostatic.advective_cfl");
  const bool has_semi_implicit_key = std::ranges::any_of(
      seen_keys, [](const auto& key) { return key.starts_with("semi_implicit."); });
  if (config.kind != ExperimentKind::kDryHydrostatic &&
      (has_time_integrator || has_advective_cfl || has_semi_implicit_key))
    throw std::runtime_error("semi-implicit keys are valid only for dry_hydrostatic");
  if (config.kind == ExperimentKind::kDryHydrostatic) {
    if (has_time_integrator) {
      if (!has_advective_cfl)
        throw std::runtime_error(
            "semi-implicit integration requires dry_hydrostatic.advective_cfl");
      require_keys(seen_keys, kSemiImplicitRequiredKeys);
    } else if (has_advective_cfl || has_semi_implicit_key) {
      throw std::runtime_error(
          "semi-implicit options require dry_hydrostatic.time_integrator");
    }
  }

  const bool has_orography_kind = seen_keys.contains("orography.kind");
  const bool has_input_file = seen_keys.contains("orography.input_file");
  const bool has_input_fingerprint =
      seen_keys.contains("orography.input_fingerprint_fnv1a64");
  const bool has_smoothing = seen_keys.contains("orography.smoothing_passes");
  if (!has_orography_kind && (has_input_file || has_input_fingerprint || has_smoothing))
    throw std::runtime_error("orography options require orography.kind");
  if (config.orography.kind == OrographyKind::kLatLonCsv &&
      !(has_input_file && has_input_fingerprint && has_smoothing))
    throw std::runtime_error(
        "latlon_csv requires input, fingerprint, and smoothing keys");
  if (config.orography.kind != OrographyKind::kLatLonCsv &&
      (has_input_file || has_input_fingerprint || has_smoothing))
    throw std::runtime_error("orography input options are valid only for latlon_csv");

  const bool has_surface_geography = seen_keys.contains("surface.geography");
  const bool has_uniform_fraction = seen_keys.contains("surface.uniform_land_fraction");
  const bool has_surface_input = seen_keys.contains("surface.input_file");
  const bool has_surface_fingerprint =
      seen_keys.contains("surface.input_fingerprint_fnv1a64");
  const bool has_surface_quadrature = seen_keys.contains("surface.quadrature_order");
  const bool has_surface_smoothing = seen_keys.contains("surface.smoothing_passes");
  if (!has_surface_geography &&
      (has_uniform_fraction || has_surface_input || has_surface_fingerprint ||
       has_surface_quadrature || has_surface_smoothing))
    throw std::runtime_error("surface options require surface.geography");
  if (has_surface_geography &&
      config.surface->geography == SurfaceGeography::kUniform && !has_uniform_fraction)
    throw std::runtime_error(
        "uniform surface geography requires surface.uniform_land_fraction");
  if (has_surface_geography && config.surface->geography == SurfaceGeography::kEarth &&
      has_uniform_fraction)
    throw std::runtime_error(
        "earth surface geography rejects surface.uniform_land_fraction");
  if (has_surface_geography && config.surface->geography == SurfaceGeography::kEarth &&
      !(has_surface_input && has_surface_fingerprint && has_surface_quadrature &&
        has_surface_smoothing))
    throw std::runtime_error(
        "earth surface geography requires input, fingerprint, quadrature, and "
        "smoothing keys");
  if (has_surface_geography &&
      config.surface->geography == SurfaceGeography::kUniform &&
      (has_surface_input || has_surface_fingerprint || has_surface_quadrature ||
       has_surface_smoothing))
    throw std::runtime_error("uniform surface geography rejects Earth data options");

  const bool has_forcing_geometry = seen_keys.contains("forcing.geometry");
  if ((config.physics.kind == PhysicsKind::kPlanetaryNewtonian) != has_forcing_geometry)
    throw std::runtime_error(
        "planetary_newtonian physics and forcing.geometry must be specified together");
  if (config.physics.kind == PhysicsKind::kPlanetaryNewtonian &&
      config.physics.geometry == ForcingGeometry::kSubstellar && !has_orbit_key)
    throw std::runtime_error("substellar forcing requires all orbit keys");
  const bool has_surface_physics_key =
      std::ranges::any_of(kSurfacePhysicsRequiredKeys,
                          [&](const auto key) { return seen_keys.contains(key); });
  if (config.physics.kind == PhysicsKind::kSurfaceEnergyBalance ||
      config.physics.kind == PhysicsKind::kGrayRadiation) {
    require_keys(seen_keys, kSurfacePhysicsRequiredKeys);
    if (!has_orbit_key || !has_surface_geography)
      throw std::runtime_error(
          "surface radiation physics requires orbit and surface geography keys");
  } else if (has_surface_physics_key) {
    throw std::runtime_error(
        "surface physics parameters require surface radiation physics");
  }
  const bool has_radiation_key = std::ranges::any_of(
      kRadiationRequiredKeys, [&](const auto key) { return seen_keys.contains(key); });
  if (config.physics.kind == PhysicsKind::kGrayRadiation) {
    require_keys(seen_keys, kRadiationRequiredKeys);
  } else if (has_radiation_key) {
    throw std::runtime_error("radiation parameters require gray_radiation");
  }

  config.validate();
  return config;
}

ExperimentConfig load_experiment_config(const std::filesystem::path& path) {
  std::ifstream input(path);
  input.imbue(std::locale::classic());
  if (!input) {
    throw std::runtime_error("unable to open configuration file: " + path.string());
  }
  auto config = parse_experiment_config(input);
  config.source_directory = path.parent_path();
  return config;
}

void write_experiment_config(std::ostream& output, const ExperimentConfig& config) {
  config.validate();
  output.imbue(std::locale::classic());
  output << std::setprecision(std::numeric_limits<Real>::max_digits10)
         << "experiment.kind = " << experiment_kind_name(config.kind) << '\n'
         << "planet.radius_m = " << config.planet.radius_m << '\n'
         << "planet.rotation_rate_rad_s = " << config.planet.rotation_rate_rad_s << '\n'
         << "planet.gravity_m_s2 = " << config.planet.gravity_m_s2 << '\n'
         << "planet.gas_constant_j_kg_k = " << config.planet.gas_constant_j_kg_k << '\n'
         << "planet.heat_capacity_cp_j_kg_k = " << config.planet.heat_capacity_cp_j_kg_k
         << '\n'
         << "planet.reference_pressure_pa = " << config.planet.reference_pressure_pa
         << '\n'
         << "run.start_time_s = " << config.run.start_time_s << '\n'
         << "run.end_time_s = " << config.run.end_time_s << '\n'
         << "run.time_step_s = " << config.run.time_step_s << '\n'
         << "run.random_seed = " << config.run.random_seed << '\n';
  if (config.orbit.has_value())
    output << "orbit.period_s = " << config.orbit->period_s << '\n'
           << "orbit.eccentricity = " << config.orbit->eccentricity << '\n'
           << "orbit.obliquity_rad = " << config.orbit->obliquity_rad << '\n'
           << "orbit.longitude_of_periapsis_rad = "
           << config.orbit->longitude_of_periapsis_rad << '\n'
           << "orbit.initial_mean_anomaly_rad = "
           << config.orbit->initial_mean_anomaly_rad << '\n'
           << "orbit.initial_substellar_longitude_rad = "
           << config.orbit->initial_substellar_longitude_rad << '\n'
           << "star.flux_at_semimajor_axis_w_m2 = "
           << config.orbit->stellar_flux_at_semimajor_axis_w_m2 << '\n';
  if (config.kind == ExperimentKind::kOde) {
    output << "ode.initial_value = " << config.ode.initial_value << '\n'
           << "ode.decay_rate_s_1 = " << config.ode.decay_rate_s_1 << '\n';
  } else if (config.kind == ExperimentKind::kSphereTransport) {
    output << "grid.cells_per_panel = " << config.grid.cells_per_panel << '\n'
           << "transport.test_case = "
           << transport_test_case_name(config.transport.test_case) << '\n'
           << "transport.initial_condition = "
           << initial_condition_name(config.transport.initial_condition) << '\n'
           << "transport.scheme = " << transport_scheme_name(config.transport.scheme)
           << '\n'
           << "transport.limiter = " << limiter_name(config.transport.limiter) << '\n'
           << "transport.cfl = " << config.transport.cfl << '\n'
           << "transport.rotation_axis_x = " << config.transport.rotation_axis_x << '\n'
           << "transport.rotation_axis_y = " << config.transport.rotation_axis_y << '\n'
           << "transport.rotation_axis_z = " << config.transport.rotation_axis_z << '\n'
           << "transport.angular_speed_rad_s = " << config.transport.angular_speed_rad_s
           << '\n';
  } else if (config.kind == ExperimentKind::kShallowWater) {
    output << "grid.cells_per_panel = " << config.grid.cells_per_panel << '\n'
           << "shallow_water.test_case = "
           << shallow_water_test_case_name(config.shallow_water.test_case) << '\n'
           << "shallow_water.scheme = "
           << shallow_water_scheme_name(config.shallow_water.scheme) << '\n'
           << "shallow_water.reconstruction = "
           << reconstruction_name(config.shallow_water.reconstruction) << '\n'
           << "shallow_water.limiter = " << limiter_name(config.shallow_water.limiter)
           << '\n'
           << "shallow_water.cfl = " << config.shallow_water.cfl << '\n'
           << "shallow_water.mean_depth_m = " << config.shallow_water.mean_depth_m
           << '\n'
           << "shallow_water.depth_floor_m = " << config.shallow_water.depth_floor_m
           << '\n'
           << "shallow_water.diffusion_kind = "
           << diffusion_kind_name(config.shallow_water.diffusion_kind) << '\n'
           << "shallow_water.diffusion_coefficient = "
           << config.shallow_water.diffusion_coefficient << '\n'
           << "shallow_water.flow_axis_x = " << config.shallow_water.flow_axis_x << '\n'
           << "shallow_water.flow_axis_y = " << config.shallow_water.flow_axis_y << '\n'
           << "shallow_water.flow_axis_z = " << config.shallow_water.flow_axis_z << '\n'
           << "shallow_water.maximum_velocity_m_s = "
           << config.shallow_water.maximum_velocity_m_s << '\n'
           << "diagnostics.interval_steps = " << config.diagnostics.interval_steps
           << '\n';
  } else {
    const auto write_list = [&output](const std::string_view key,
                                      const std::vector<Real>& values) {
      output << key << " = ";
      for (std::size_t index = 0; index < values.size(); ++index) {
        if (index != 0) {
          output << ',';
        }
        output << values[index];
      }
      output << '\n';
    };
    if (config.kind == ExperimentKind::kVerticalColumn) {
      output << "vertical.test_case = "
             << vertical_test_case_name(config.vertical.test_case) << '\n';
    }
    output << "vertical.levels = " << config.vertical.levels << '\n';
    write_list("vertical.a_half_pa", config.vertical.a_half_pa);
    write_list("vertical.b_half", config.vertical.b_half);
    output << "vertical.surface_pressure_pa = " << config.vertical.surface_pressure_pa
           << '\n'
           << "vertical.minimum_surface_pressure_pa = "
           << config.vertical.minimum_surface_pressure_pa << '\n'
           << "vertical.maximum_surface_pressure_pa = "
           << config.vertical.maximum_surface_pressure_pa << '\n'
           << "vertical.minimum_pressure_thickness_pa = "
           << config.vertical.minimum_pressure_thickness_pa << '\n'
           << "vertical.surface_geopotential_m2_s2 = "
           << config.vertical.surface_geopotential_m2_s2 << '\n'
           << "vertical.initial_temperature_k = "
           << config.vertical.initial_temperature_k << '\n'
           << "vertical.initial_potential_temperature_k = "
           << config.vertical.initial_potential_temperature_k << '\n'
           << "vertical.temperature_floor_k = " << config.vertical.temperature_floor_k
           << '\n'
           << "vertical.transport_scheme = "
           << vertical_transport_scheme_name(config.vertical.transport_scheme) << '\n'
           << "vertical.limiter = " << vertical_limiter_name(config.vertical.limiter)
           << '\n'
           << "vertical.cfl = " << config.vertical.cfl << '\n'
           << "diagnostics.interval_steps = " << config.diagnostics.interval_steps
           << '\n';
    if (config.kind == ExperimentKind::kVerticalColumn) {
      output << "vertical.forcing_amplitude = " << config.vertical.forcing_amplitude
             << '\n';
    }
    if (config.kind == ExperimentKind::kDryHydrostatic) {
      output << "grid.cells_per_panel = " << config.grid.cells_per_panel << '\n'
             << "dry_hydrostatic.test_case = "
             << dry_hydrostatic_test_case_name(config.dry_hydrostatic.test_case) << '\n'
             << "dry_hydrostatic.reconstruction = "
             << reconstruction_name(config.dry_hydrostatic.reconstruction) << '\n'
             << "dry_hydrostatic.limiter = "
             << limiter_name(config.dry_hydrostatic.limiter) << '\n'
             << "dry_hydrostatic.cfl = " << config.dry_hydrostatic.cfl << '\n'
             << "dry_hydrostatic.diffusion_kind = "
             << diffusion_kind_name(config.dry_hydrostatic.diffusion_kind) << '\n'
             << "dry_hydrostatic.diffusion_coefficient = "
             << config.dry_hydrostatic.diffusion_coefficient << '\n';
      if (config.dry_hydrostatic.time_integrator ==
          DryHydrostaticTimeIntegrator::kSemiImplicit) {
        const auto& semi_implicit = *config.semi_implicit;
        output << "dry_hydrostatic.time_integrator = "
               << dry_hydrostatic_time_integrator_name(
                      config.dry_hydrostatic.time_integrator)
               << '\n'
               << "dry_hydrostatic.advective_cfl = "
               << config.dry_hydrostatic.advective_cfl << '\n'
               << "semi_implicit.reference_surface_pressure_pa = "
               << semi_implicit.reference_surface_pressure_pa << '\n'
               << "semi_implicit.reference_temperature_k = "
               << semi_implicit.reference_temperature_k << '\n'
               << "semi_implicit.reference_update = "
               << (semi_implicit.reference_update ==
                           SemiImplicitReferenceUpdate::kPerStep
                       ? "per_step"
                       : "fixed")
               << '\n'
               << "semi_implicit.implicit_weight = " << semi_implicit.implicit_weight
               << '\n'
               << "semi_implicit.wave_cfl_threshold = "
               << semi_implicit.wave_cfl_threshold << '\n'
               << "semi_implicit.maximum_implicit_modes = "
               << semi_implicit.maximum_implicit_modes << '\n'
               << "semi_implicit.nonlinear_iterations = "
               << semi_implicit.nonlinear_iterations << '\n'
               << "semi_implicit.linear_relative_tolerance = "
               << semi_implicit.linear_relative_tolerance << '\n'
               << "semi_implicit.linear_absolute_tolerance = "
               << semi_implicit.linear_absolute_tolerance << '\n'
               << "semi_implicit.linear_maximum_iterations = "
               << semi_implicit.linear_maximum_iterations << '\n'
               << "semi_implicit.gmres_restart = " << semi_implicit.gmres_restart
               << '\n'
               << "semi_implicit.minimum_time_step_s = "
               << semi_implicit.minimum_time_step_s << '\n';
      }
      if (config.physics.kind != PhysicsKind::kNone)
        output << "physics.kind = " << physics_kind_name(config.physics.kind) << '\n';
      if (config.physics.kind == PhysicsKind::kPlanetaryNewtonian)
        output << "forcing.geometry = "
               << forcing_geometry_name(config.physics.geometry) << '\n';
      if (config.surface.has_value()) {
        output << "surface.geography = "
               << surface_geography_name(config.surface->geography) << '\n';
        if (config.surface->geography == SurfaceGeography::kUniform)
          output << "surface.uniform_land_fraction = "
                 << config.surface->uniform_land_fraction << '\n';
        else
          output << "surface.input_file = " << config.surface->input_file << '\n'
                 << "surface.input_fingerprint_fnv1a64 = "
                 << config.surface->input_fingerprint_fnv1a64 << '\n'
                 << "surface.quadrature_order = " << config.surface->quadrature_order
                 << '\n'
                 << "surface.smoothing_passes = " << config.surface->smoothing_passes
                 << '\n';
        if (config.physics.kind == PhysicsKind::kSurfaceEnergyBalance ||
            config.physics.kind == PhysicsKind::kGrayRadiation)
          output << "surface.land_heat_capacity_j_m2_k = "
                 << config.surface->land_heat_capacity_j_m2_k << '\n'
                 << "surface.ocean_heat_capacity_j_m2_k = "
                 << config.surface->ocean_heat_capacity_j_m2_k << '\n'
                 << "surface.initial_temperature_k = "
                 << config.surface->initial_temperature_k << '\n'
                 << "surface.albedo = " << config.surface->albedo << '\n'
                 << "surface.emissivity = " << config.surface->emissivity << '\n'
                 << "surface.air_exchange_coefficient_w_m2_k = "
                 << config.surface->air_exchange_coefficient_w_m2_k << '\n'
                 << "surface.internal_heat_flux_w_m2 = "
                 << config.surface->internal_heat_flux_w_m2 << '\n'
                 << "surface.cfl = " << config.surface->cfl << '\n';
      }
      if (config.radiation.has_value())
        output << "radiation.shortwave_absorption_m2_kg = "
               << config.radiation->shortwave_absorption_m2_kg << '\n'
               << "radiation.longwave_absorption_ref_m2_kg = "
               << config.radiation->longwave_absorption_ref_m2_kg << '\n'
               << "radiation.reference_pressure_pa = "
               << config.radiation->reference_pressure_pa << '\n'
               << "radiation.longwave_pressure_exponent = "
               << config.radiation->longwave_pressure_exponent << '\n'
               << "radiation.longwave_diffusivity_factor = "
               << config.radiation->longwave_diffusivity_factor << '\n'
               << "radiation.shortwave_diffuse_factor = "
               << config.radiation->shortwave_diffuse_factor << '\n'
               << "radiation.cfl = " << config.radiation->cfl << '\n';
    }
    if (config.orography.kind != OrographyKind::kFlat) {
      output << "orography.kind = " << orography_kind_name(config.orography.kind)
             << '\n';
      if (config.orography.kind == OrographyKind::kLatLonCsv)
        output << "orography.input_file = " << config.orography.input_file << '\n'
               << "orography.input_fingerprint_fnv1a64 = "
               << config.orography.input_fingerprint_fnv1a64 << '\n'
               << "orography.smoothing_passes = " << config.orography.smoothing_passes
               << '\n';
    }
  }
  output << "output.directory = " << config.output_directory << '\n';

  if (!output) {
    throw std::runtime_error("failed while writing configuration stream");
  }
}

std::string_view experiment_kind_name(const ExperimentKind kind) noexcept {
  switch (kind) {
    case ExperimentKind::kOde:
      return "ode";
    case ExperimentKind::kSphereTransport:
      return "sphere_transport";
    case ExperimentKind::kShallowWater:
      return "shallow_water";
    case ExperimentKind::kVerticalColumn:
      return "vertical_column";
    case ExperimentKind::kDryHydrostatic:
      return "dry_hydrostatic";
  }
  return "unknown";
}

std::string_view dry_hydrostatic_test_case_name(
    const DryHydrostaticTestCase value) noexcept {
  switch (value) {
    case DryHydrostaticTestCase::kIsothermalRest:
      return "isothermal_rest";
    case DryHydrostaticTestCase::kSolidBodyTransport:
      return "solid_body_transport";
    case DryHydrostaticTestCase::kDcmipDeformational:
      return "dcmip_deformational";
    case DryHydrostaticTestCase::kDcmipHadley:
      return "dcmip_hadley";
    case DryHydrostaticTestCase::kLinearWave:
      return "linear_wave";
    case DryHydrostaticTestCase::kDcmip200Rest:
      return "dcmip_2_0_0_rest";
    case DryHydrostaticTestCase::kLinearMountainWave:
      return "linear_mountain_wave";
    case DryHydrostaticTestCase::kJw06Steady:
      return "jw06_steady";
    case DryHydrostaticTestCase::kJw06Baroclinic:
      return "jw06_baroclinic";
    case DryHydrostaticTestCase::kUmjs14Steady:
      return "umjs14_steady";
    case DryHydrostaticTestCase::kUmjs14Baroclinic:
      return "umjs14_baroclinic";
    case DryHydrostaticTestCase::kHeldSuarez:
      return "held_suarez";
  }
  return "unknown";
}

std::string_view dry_hydrostatic_time_integrator_name(
    const DryHydrostaticTimeIntegrator integrator) noexcept {
  switch (integrator) {
    case DryHydrostaticTimeIntegrator::kExplicitSspRk3:
      return "explicit_ssprk3";
    case DryHydrostaticTimeIntegrator::kSemiImplicit:
      return "semi_implicit";
  }
  return "unknown";
}

std::string_view orography_kind_name(const OrographyKind kind) noexcept {
  switch (kind) {
    case OrographyKind::kFlat:
      return "flat";
    case OrographyKind::kDcmip200:
      return "dcmip_2_0_0";
    case OrographyKind::kWilliamson5:
      return "williamson5";
    case OrographyKind::kLinearBell:
      return "linear_bell";
    case OrographyKind::kJw06:
      return "jw06";
    case OrographyKind::kLatLonCsv:
      return "latlon_csv";
  }
  return "unknown";
}

std::string_view physics_kind_name(const PhysicsKind kind) noexcept {
  switch (kind) {
    case PhysicsKind::kNone:
      return "none";
    case PhysicsKind::kHeldSuarez:
      return "held_suarez";
    case PhysicsKind::kPlanetaryNewtonian:
      return "planetary_newtonian";
    case PhysicsKind::kSurfaceEnergyBalance:
      return "surface_energy_balance";
    case PhysicsKind::kGrayRadiation:
      return "gray_radiation";
  }
  return "unknown";
}

std::string_view forcing_geometry_name(const ForcingGeometry geometry) noexcept {
  return geometry == ForcingGeometry::kAxisymmetric ? "axisymmetric" : "substellar";
}

std::string_view surface_geography_name(const SurfaceGeography geography) noexcept {
  switch (geography) {
    case SurfaceGeography::kUniform:
      return "uniform";
    case SurfaceGeography::kEarth:
      return "earth";
  }
  return "unknown";
}

std::string_view vertical_test_case_name(const VerticalTestCase test_case) noexcept {
  switch (test_case) {
    case VerticalTestCase::kIsothermal:
      return "isothermal";
    case VerticalTestCase::kDryAdiabatic:
      return "dry_adiabatic";
    case VerticalTestCase::kMovingSurfacePressure:
      return "moving_surface_pressure";
    case VerticalTestCase::kManufacturedTransport:
      return "manufactured_transport";
  }
  return "unknown";
}

std::string_view vertical_transport_scheme_name(
    const VerticalTransportScheme scheme) noexcept {
  return scheme == VerticalTransportScheme::kDonorCell ? "donor_cell" : "linear";
}

std::string_view vertical_limiter_name(const VerticalLimiterKind limiter) noexcept {
  return limiter == VerticalLimiterKind::kNone ? "none" : "minmod";
}

std::string_view shallow_water_scheme_name(const ShallowWaterScheme scheme) noexcept {
  return scheme == ShallowWaterScheme::kRusanov ? "rusanov" : "compatible";
}

std::string_view reconstruction_name(const ReconstructionKind reconstruction) noexcept {
  return reconstruction == ReconstructionKind::kPiecewiseConstant ? "piecewise_constant"
                                                                  : "linear";
}

std::string_view shallow_water_test_case_name(
    const ShallowWaterTestCase test_case) noexcept {
  switch (test_case) {
    case ShallowWaterTestCase::kRest:
      return "rest";
    case ShallowWaterTestCase::kLinearWave:
      return "linear_wave";
    case ShallowWaterTestCase::kGeostrophicAdjustment:
      return "geostrophic_adjustment";
    case ShallowWaterTestCase::kWilliamson2:
      return "williamson2";
    case ShallowWaterTestCase::kWilliamson5:
      return "williamson5";
    case ShallowWaterTestCase::kWilliamson6:
      return "williamson6";
    case ShallowWaterTestCase::kGalewsky:
      return "galewsky";
  }
  return "unknown";
}

std::string_view diffusion_kind_name(const DiffusionKind kind) noexcept {
  switch (kind) {
    case DiffusionKind::kNone:
      return "none";
    case DiffusionKind::kLaplacian:
      return "laplacian";
    case DiffusionKind::kBiharmonic:
      return "biharmonic";
  }
  return "unknown";
}

std::string_view transport_scheme_name(const TransportScheme scheme) noexcept {
  return scheme == TransportScheme::kUpwind ? "upwind" : "linear";
}

std::string_view limiter_name(const LimiterKind limiter) noexcept {
  return limiter == LimiterKind::kNone ? "none" : "barth_jespersen";
}

std::string_view transport_test_case_name(const TransportTestCase test_case) noexcept {
  switch (test_case) {
    case TransportTestCase::kSolidBody:
      return "solid_body";
    case TransportTestCase::kDeformational:
      return "deformational";
    case TransportTestCase::kDivergent:
      return "divergent";
  }
  return "unknown";
}

std::string_view initial_condition_name(
    const InitialConditionKind initial_condition) noexcept {
  switch (initial_condition) {
    case InitialConditionKind::kConstant:
      return "constant";
    case InitialConditionKind::kGaussianHill:
      return "gaussian_hill";
    case InitialConditionKind::kCosineBell:
      return "cosine_bell";
    case InitialConditionKind::kSlottedCylinder:
      return "slotted_cylinder";
  }
  return "unknown";
}

}  // namespace mps
