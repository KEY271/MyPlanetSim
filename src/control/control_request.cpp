#include "myplanetsim/control/control_request.hpp"

#include <array>
#include <cctype>
#include <charconv>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <locale>
#include <set>
#include <stdexcept>
#include <string_view>
#include <vector>

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
  } else if (key == "control.cells_per_panel") {
    request.cells_per_panel = parse_uint(value, line, key);
  } else if (key == "control.end_time_s") {
    request.end_time_s = parse_real(value, line, key);
  } else if (key == "control.maximum_time_step_s") {
    request.maximum_time_step_s = parse_real(value, line, key);
  } else if (key == "control.frame_interval_steps") {
    request.frame_interval_steps = parse_uint(value, line, key);
  } else if (key == "control.frame_directory") {
    request.frame_directory = std::filesystem::path(value);
  } else if (key == "initial_edits.count") {
    const auto count = parse_uint(value, line, key);
    if (count > kControlMaxEditCount) {
      throw parse_error(line, "initial_edits.count exceeds the supported limit");
    }
    request.initial_edits.resize(static_cast<std::size_t>(count));
  } else if (key.starts_with("initial_edits.")) {
    const auto first_dot = key.find('.', std::string_view("initial_edits.").size());
    if (first_dot == std::string_view::npos) {
      throw parse_error(line, "invalid initial edit key");
    }
    const auto index_text =
        key.substr(std::string_view("initial_edits.").size(),
                   first_dot - std::string_view("initial_edits.").size());
    const auto index = parse_uint(index_text, line, key);
    if (index >= kControlMaxEditCount) {
      throw parse_error(line, "initial edit index exceeds the supported limit");
    }
    if (index >= request.initial_edits.size()) {
      request.initial_edits.resize(static_cast<std::size_t>(index + 1));
    }
    auto& edit = request.initial_edits[static_cast<std::size_t>(index)];
    const auto field = key.substr(first_dot + 1);
    if (field == "kind") {
      if (value != "gaussian_depth") {
        throw parse_error(line, "unknown initial edit kind " + std::string(value));
      }
      edit.kind = InitialConditionEditKind::kGaussianDepth;
    } else if (field == "center_x") {
      edit.center_unit.x = parse_real(value, line, key);
    } else if (field == "center_y") {
      edit.center_unit.y = parse_real(value, line, key);
    } else if (field == "center_z") {
      edit.center_unit.z = parse_real(value, line, key);
    } else if (field == "amplitude_m") {
      edit.amplitude_m = parse_real(value, line, key);
    } else if (field == "sigma_rad") {
      edit.sigma_rad = parse_real(value, line, key);
    } else if (field == "mass_policy") {
      if (value == "preserve_global") {
        edit.mass_policy = MassPolicy::kPreserveGlobal;
      } else if (value == "allow_change") {
        edit.mass_policy = MassPolicy::kAllowChange;
      } else {
        throw parse_error(line,
                          "unknown initial edit mass policy " + std::string(value));
      }
    } else {
      throw parse_error(line, "unknown initial edit field " + std::string(field));
    }
  } else {
    throw parse_error(line, "unknown key " + std::string(key));
  }
}

[[nodiscard]] bool is_run_id_character(const char character) {
  return std::isalnum(static_cast<unsigned char>(character)) != 0 || character == '_' ||
         character == '-';
}

[[nodiscard]] bool has_path_escape(const std::filesystem::path& path) {
  if (path.empty() || path.is_absolute() || path.has_root_name() ||
      path.has_root_directory()) {
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
  if (request.cells_per_panel == 0 || request.cells_per_panel > kControlMaxCellsPerPanel) {
    throw std::invalid_argument("control.cells_per_panel must be in [1, 96]");
  }
  require_finite(request.end_time_s, "control.end_time_s");
  if (!(request.end_time_s > 0.0 && request.end_time_s <= 31536000.0)) {
    throw std::invalid_argument("control.end_time_s must be in (0, 31536000]");
  }
  require_finite(request.maximum_time_step_s, "control.maximum_time_step_s");
  if (!(request.maximum_time_step_s > 0.0 && request.maximum_time_step_s <= 86400.0)) {
    throw std::invalid_argument("control.maximum_time_step_s must be in (0, 86400]");
  }
  if (request.frame_interval_steps == 0 || request.frame_interval_steps > 1000000) {
    throw std::invalid_argument("control.frame_interval_steps must be in [1, 1000000]");
  }
  // Every published frame is a file on disk and a decoded dataset in the client, so the
  // request is rejected when even the largest permitted time step would exceed the budget.
  const Real estimated_frames =
      std::floor(std::ceil(request.end_time_s / request.maximum_time_step_s) /
                 static_cast<Real>(request.frame_interval_steps)) +
      2.0;
  if (estimated_frames > static_cast<Real>(kControlMaxPublishedFrames)) {
    throw std::invalid_argument(
        "control request would publish more than the supported number of frames");
  }
  if (has_path_escape(request.frame_directory)) {
    throw std::invalid_argument("control.frame_directory must be a safe relative path");
  }
  if (request.initial_edits.size() > kControlMaxEditCount) {
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
  static constexpr std::array<std::string_view, 8> required_keys{
      "control.format_version",      "control.run_id",
      "control.cells_per_panel",     "control.end_time_s",
      "control.maximum_time_step_s", "control.frame_interval_steps",
      "control.frame_directory",     "initial_edits.count"};
  for (const auto key : required_keys) {
    if (!seen_keys.contains(key)) {
      throw std::runtime_error("control request is missing required key " +
                               std::string(key));
    }
  }
  for (std::size_t index = 0; index < request.initial_edits.size(); ++index) {
    const std::string prefix = "initial_edits." + std::to_string(index) + ".";
    static constexpr std::array<std::string_view, 7> edit_fields{
        "kind",        "center_x",  "center_y",   "center_z",
        "amplitude_m", "sigma_rad", "mass_policy"};
    for (const auto field : edit_fields) {
      if (!seen_keys.contains(prefix + std::string(field))) {
        throw std::runtime_error("control request is missing required key " + prefix +
                                 std::string(field));
      }
    }
  }
  for (const auto& key : seen_keys) {
    if (!key.starts_with("initial_edits.") || key == "initial_edits.count") {
      continue;
    }
    const auto first_dot = key.find('.', std::string_view("initial_edits.").size());
    const auto index_text =
        key.substr(std::string_view("initial_edits.").size(),
                   first_dot - std::string_view("initial_edits.").size());
    std::uint64_t index = 0;
    const auto result = std::from_chars(index_text.data(),
                                        index_text.data() + index_text.size(), index);
    if (result.ec != std::errc{} ||
        result.ptr != index_text.data() + index_text.size() ||
        index >= request.initial_edits.size()) {
      throw std::runtime_error("initial edit key is outside initial_edits.count");
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
         << "control.cells_per_panel = " << request.cells_per_panel << '\n'
         << "control.end_time_s = " << request.end_time_s << '\n'
         << "control.maximum_time_step_s = " << request.maximum_time_step_s << '\n'
         << "control.frame_interval_steps = " << request.frame_interval_steps << '\n'
         << "control.frame_directory = " << request.frame_directory.generic_string()
         << '\n'
         << "initial_edits.count = " << request.initial_edits.size() << '\n';
  for (std::size_t index = 0; index < request.initial_edits.size(); ++index) {
    const auto& edit = request.initial_edits[index];
    output << "initial_edits." << index
           << ".kind = " << initial_condition_edit_kind_name(edit.kind) << '\n'
           << "initial_edits." << index << ".center_x = " << edit.center_unit.x << '\n'
           << "initial_edits." << index << ".center_y = " << edit.center_unit.y << '\n'
           << "initial_edits." << index << ".center_z = " << edit.center_unit.z << '\n'
           << "initial_edits." << index << ".amplitude_m = " << edit.amplitude_m << '\n'
           << "initial_edits." << index << ".sigma_rad = " << edit.sigma_rad << '\n'
           << "initial_edits." << index
           << ".mass_policy = " << mass_policy_name(edit.mass_policy) << '\n';
  }
  if (!output) {
    throw std::runtime_error("failed while writing control request");
  }
}

ExperimentConfig apply_control_request(const ExperimentConfig& config,
                                       const ControlRequestV1& request) {
  validate_control_request(request);
  ExperimentConfig resolved = config;
  resolved.grid.cells_per_panel = static_cast<Index>(request.cells_per_panel);
  resolved.run.end_time_s = request.end_time_s;
  resolved.run.time_step_s = request.maximum_time_step_s;
  resolved.diagnostics.interval_steps = request.frame_interval_steps;
  resolved.output_directory = request.frame_directory.generic_string();
  resolved.validate();
  return resolved;
}

std::string_view initial_condition_edit_kind_name(
    const InitialConditionEditKind kind) noexcept {
  switch (kind) {
    case InitialConditionEditKind::kGaussianDepth:
      return "gaussian_depth";
  }
  return "unknown";
}

std::string_view mass_policy_name(const MassPolicy policy) noexcept {
  switch (policy) {
    case MassPolicy::kPreserveGlobal:
      return "preserve_global";
    case MassPolicy::kAllowChange:
      return "allow_change";
  }
  return "unknown";
}

}  // namespace mps
