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

constexpr std::array<std::string_view, 12> kTransportRequiredKeys{
    "experiment.kind",
    "grid.cells_per_panel",
    "grid.halo_width",
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

void assign_value(ExperimentConfig& config, const std::string_view key,
                  const std::string_view value, const std::size_t line) {
  if (key == "experiment.kind") {
    if (value == "ode") {
      config.kind = ExperimentKind::kOde;
    } else if (value == "sphere_transport") {
      config.kind = ExperimentKind::kSphereTransport;
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
  } else if (key == "grid.halo_width") {
    config.grid.halo_width = parse_index(value, line, key);
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
  } else {
    if (grid.cells_per_panel <= 0) {
      throw std::invalid_argument("grid.cells_per_panel must be positive");
    }
    if (grid.halo_width < 1) {
      throw std::invalid_argument("grid.halo_width must be at least one");
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
      if (key.starts_with("grid.") || key.starts_with("transport.")) {
        throw std::runtime_error("key " + key + " is not valid for ode experiment");
      }
    }
  } else {
    require_keys(seen_keys, kTransportRequiredKeys);
    for (const auto& key : seen_keys) {
      if (key.starts_with("ode.")) {
        throw std::runtime_error("key " + key +
                                 " is not valid for sphere_transport experiment");
      }
    }
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
  } else {
    output << "grid.cells_per_panel = " << config.grid.cells_per_panel << '\n'
           << "grid.halo_width = " << config.grid.halo_width << '\n'
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
  }
  output << "output.directory = " << config.output_directory << '\n';

  if (!output) {
    throw std::runtime_error("failed while writing configuration stream");
  }
}

std::string_view experiment_kind_name(const ExperimentKind kind) noexcept {
  return kind == ExperimentKind::kOde ? "ode" : "sphere_transport";
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
