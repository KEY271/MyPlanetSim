#include "myplanetsim/config/experiment_config.hpp"

#include <array>
#include <cctype>
#include <charconv>
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

constexpr std::array<std::string_view, 13> kRequiredKeys{
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
    "ode.initial_value",
    "ode.decay_rate_s_1",
    "output.directory",
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

void assign_value(ExperimentConfig& config, const std::string_view key,
                  const std::string_view value, const std::size_t line) {
  if (key == "planet.radius_m") {
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
  } else if (key == "output.directory") {
    config.output_directory = value;
  } else {
    throw parse_error(line, "unknown key " + std::string(key));
  }
}

}  // namespace

void ExperimentConfig::validate() const {
  planet.validate();
  require_non_negative(run.start_time_s, "run.start_time_s");
  require_finite(run.end_time_s, "run.end_time_s");
  require_positive(run.time_step_s, "run.time_step_s");
  require_finite(ode.initial_value, "ode.initial_value");
  require_non_negative(ode.decay_rate_s_1, "ode.decay_rate_s_1");

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

  for (const auto key : kRequiredKeys) {
    if (!seen_keys.contains(key)) {
      throw std::runtime_error("configuration is missing required key " +
                               std::string(key));
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
         << "run.random_seed = " << config.run.random_seed << '\n'
         << "ode.initial_value = " << config.ode.initial_value << '\n'
         << "ode.decay_rate_s_1 = " << config.ode.decay_rate_s_1 << '\n'
         << "output.directory = " << config.output_directory << '\n';

  if (!output) {
    throw std::runtime_error("failed while writing configuration stream");
  }
}

}  // namespace mps
