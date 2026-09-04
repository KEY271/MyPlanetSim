#include "myplanetsim/config/experiment_config.hpp"

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
      config.dry_hydrostatic.test_case =
          DryHydrostaticTestCase::kLinearMountainWave;
    else if (value == "umjs14_steady")
      config.dry_hydrostatic.test_case = DryHydrostaticTestCase::kUmjs14Steady;
    else if (value == "umjs14_baroclinic")
      config.dry_hydrostatic.test_case = DryHydrostaticTestCase::kUmjs14Baroclinic;
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
    if ((dry_hydrostatic.test_case ==
         DryHydrostaticTestCase::kLinearMountainWave) !=
        (orography.kind == OrographyKind::kLinearBell))
      throw std::invalid_argument(
          "linear_mountain_wave requires linear_bell orography and vice versa");
  }

  if (orography.kind == OrographyKind::kFlat) {
    if (!orography.input_file.empty() ||
        !orography.input_fingerprint_fnv1a64.empty() ||
        orography.smoothing_passes != 0)
      throw std::invalid_argument("flat orography does not accept input options");
  } else if (orography.kind == OrographyKind::kLatLonCsv) {
    if (orography.input_file.empty() ||
        orography.input_file.find_first_of("\r\n") != std::string::npos)
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
  if (config.kind == ExperimentKind::kOde) {
    require_keys(seen_keys, kOdeRequiredKeys);
    for (const auto& key : seen_keys) {
      if (key.starts_with("grid.") || key.starts_with("transport.") ||
          key.starts_with("shallow_water.") || key.starts_with("diagnostics.") ||
          key.starts_with("orography.")) {
        throw std::runtime_error("key " + key + " is not valid for ode experiment");
      }
    }
  } else if (config.kind == ExperimentKind::kSphereTransport) {
    require_keys(seen_keys, kTransportRequiredKeys);
    for (const auto& key : seen_keys) {
      if (key.starts_with("ode.") || key.starts_with("shallow_water.") ||
          key.starts_with("diagnostics.") || key.starts_with("orography.")) {
        throw std::runtime_error("key " + key +
                                 " is not valid for sphere_transport experiment");
      }
    }
  } else if (config.kind == ExperimentKind::kShallowWater) {
    require_keys(seen_keys, kShallowWaterRequiredKeys);
    for (const auto& key : seen_keys) {
      if (key.starts_with("ode.") || key.starts_with("transport.")) {
        throw std::runtime_error("key " + key +
                                 " is not valid for shallow_water experiment");
      }
    }
  } else if (config.kind == ExperimentKind::kVerticalColumn) {
    require_keys(seen_keys, kVerticalRequiredKeys);
    for (const auto& key : seen_keys) {
      if (key.starts_with("ode.") || key.starts_with("grid.") ||
          key.starts_with("transport.") || key.starts_with("shallow_water.") ||
          key.starts_with("orography.")) {
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


  const bool has_orography_kind = seen_keys.contains("orography.kind");
  const bool has_input_file = seen_keys.contains("orography.input_file");
  const bool has_input_fingerprint =
      seen_keys.contains("orography.input_fingerprint_fnv1a64");
  const bool has_smoothing = seen_keys.contains("orography.smoothing_passes");
  if (!has_orography_kind && (has_input_file || has_input_fingerprint || has_smoothing))
    throw std::runtime_error("orography options require orography.kind");
  if (config.orography.kind == OrographyKind::kLatLonCsv &&
      !(has_input_file && has_input_fingerprint && has_smoothing))
    throw std::runtime_error("latlon_csv requires input, fingerprint, and smoothing keys");
  if (config.orography.kind != OrographyKind::kLatLonCsv &&
      (has_input_file || has_input_fingerprint || has_smoothing))
    throw std::runtime_error("orography input options are valid only for latlon_csv");

  config.validate();
  return config;
}

ExperimentConfig load_experiment_config(const std::filesystem::path& path) {
  std::ifstream input(path);
  input.imbue(std::locale::classic());
  if (!input) {
    throw std::runtime_error("unable to open configuration file: " + path.string());
  }
  return parse_experiment_config(input);
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
  }
  if (config.orography.kind != OrographyKind::kFlat) {
    output << "orography.kind = " << orography_kind_name(config.orography.kind)
           << '\n';
    if (config.orography.kind == OrographyKind::kLatLonCsv)
      output << "orography.input_file = " << config.orography.input_file << '\n'
             << "orography.input_fingerprint_fnv1a64 = "
             << config.orography.input_fingerprint_fnv1a64 << '\n'
             << "orography.smoothing_passes = "
             << config.orography.smoothing_passes << '\n';
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
    case DryHydrostaticTestCase::kUmjs14Steady:
      return "umjs14_steady";
    case DryHydrostaticTestCase::kUmjs14Baroclinic:
      return "umjs14_baroclinic";
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
