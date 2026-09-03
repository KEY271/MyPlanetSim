#include "myplanetsim/control/control_request.hpp"

#include <array>
#include <charconv>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <locale>
#include <set>
#include <stdexcept>
#include <string_view>

#include "myplanetsim/core/validation.hpp"

namespace mps {
namespace {

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
  return std::runtime_error("control request line " + std::to_string(line) + ": " +
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

[[nodiscard]] std::uint64_t parse_uint(const std::string_view text,
                                       const std::size_t line,
                                       const std::string_view key) {
  std::uint64_t value = 0;
  const auto result = std::from_chars(text.data(), text.data() + text.size(), value);
  if (result.ec != std::errc{} || result.ptr != text.data() + text.size()) {
    throw parse_error(line, "invalid unsigned integer value for " + std::string(key));
  }
  return value;
}

void assign_value(ControlRequestV1& request, const std::string_view key,
                  const std::string_view value, const std::size_t line) {
  if (key == "control.format_version") {
    request.format_version = parse_uint(value, line, key);
  } else if (key == "control.run_id") {
    request.run_id = std::string(value);
  } else if (key == "control.end_time_s") {
    request.end_time_s = parse_real(value, line, key);
  } else if (key == "control.maximum_time_step_s") {
    request.maximum_time_step_s = parse_real(value, line, key);
  } else if (key == "control.frame_interval_steps") {
    request.frame_interval_steps = parse_uint(value, line, key);
  } else if (key == "control.frame_directory") {
    request.frame_directory = std::filesystem::path(value);
  } else if (key == "initial_edits.count") {
    request.initial_edits_count = parse_uint(value, line, key);
  } else {
    throw parse_error(line, "unknown key " + std::string(key));
  }
}

[[nodiscard]] bool is_run_id_character(const char character) {
  return std::isalnum(static_cast<unsigned char>(character)) != 0 ||
         character == '_' || character == '-';
}

[[nodiscard]] bool has_path_escape(const std::filesystem::path& path) {
  if (path.empty() || path.is_absolute() || path.has_root_name() || path.has_root_directory()) {
    return true;
  }
  for (const auto& component : path) {
    if (component == ".." || component == ".") {
      return true;
    }
  }
  return false;
}

}  // namespace

void validate_control_request(const ControlRequestV1& request) {
  if (request.format_version != kControlProtocolVersion) {
    throw std::invalid_argument("unsupported control request format version");
  }
  if (request.run_id.empty() || request.run_id.size() > 64) {
    throw std::invalid_argument("control.run_id must contain 1 to 64 characters");
  }
  for (const char character : request.run_id) {
    if (!is_run_id_character(character)) {
      throw std::invalid_argument("control.run_id contains an invalid character");
    }
  }
  require_finite(request.end_time_s, "control.end_time_s");
  if (!(request.end_time_s > 0.0 && request.end_time_s <= 31536000.0)) {
    throw std::invalid_argument("control.end_time_s must be in (0, 31536000]");
  }
  require_finite(request.maximum_time_step_s, "control.maximum_time_step_s");
  if (!(request.maximum_time_step_s > 0.0 && request.maximum_time_step_s <= 86400.0)) {
    throw std::invalid_argument(
        "control.maximum_time_step_s must be in (0, 86400]");
  }
  if (request.frame_interval_steps == 0 || request.frame_interval_steps > 1000000) {
    throw std::invalid_argument(
        "control.frame_interval_steps must be in [1, 1000000]");
  }
  if (has_path_escape(request.frame_directory)) {
    throw std::invalid_argument("control.frame_directory must be a safe relative path");
  }
  if (request.initial_edits_count > kControlMaxEditCount) {
    throw std::invalid_argument("initial_edits.count exceeds the supported limit");
  }
}

ControlRequestV1 parse_control_request(std::istream& input) {
  ControlRequestV1 request{};
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
    assign_value(request, key, value, line_number);
  }
  if (input.bad()) {
    throw std::runtime_error("failed while reading control request stream");
  }
  static constexpr std::array<std::string_view, 7> required_keys{
      "control.format_version", "control.run_id", "control.end_time_s",
      "control.maximum_time_step_s", "control.frame_interval_steps",
      "control.frame_directory", "initial_edits.count"};
  for (const auto key : required_keys) {
    if (!seen_keys.contains(key)) {
      throw std::runtime_error("control request is missing required key " +
                               std::string(key));
    }
  }
  validate_control_request(request);
  return request;
}

ControlRequestV1 load_control_request(const std::filesystem::path& path) {
  std::ifstream input(path);
  input.imbue(std::locale::classic());
  if (!input) {
    throw std::runtime_error("unable to open control request: " + path.string());
  }
  return parse_control_request(input);
}

void write_control_request(std::ostream& output, const ControlRequestV1& request) {
  validate_control_request(request);
  output.imbue(std::locale::classic());
  output << std::setprecision(std::numeric_limits<Real>::max_digits10)
         << "control.format_version = " << request.format_version << '\n'
         << "control.run_id = " << request.run_id << '\n'
         << "control.end_time_s = " << request.end_time_s << '\n'
         << "control.maximum_time_step_s = " << request.maximum_time_step_s << '\n'
         << "control.frame_interval_steps = " << request.frame_interval_steps << '\n'
         << "control.frame_directory = " << request.frame_directory.generic_string() << '\n'
         << "initial_edits.count = " << request.initial_edits_count << '\n';
  if (!output) {
    throw std::runtime_error("failed while writing control request");
  }
}

ExperimentConfig apply_control_request(const ExperimentConfig& config,
                                       const ControlRequestV1& request) {
  validate_control_request(request);
  ExperimentConfig resolved = config;
  resolved.run.end_time_s = request.end_time_s;
  resolved.run.time_step_s = request.maximum_time_step_s;
  resolved.diagnostics.interval_steps = request.frame_interval_steps;
  resolved.output_directory = request.frame_directory.generic_string();
  resolved.validate();
  return resolved;
}

}  // namespace mps
