#include "myplanetsim/io/checkpoint.hpp"

#include <charconv>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <locale>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>

#include "myplanetsim/core/validation.hpp"

namespace mps {
namespace {

constexpr std::size_t kMaximumStateSize = 100'000'000;

[[nodiscard]] std::string read_line(std::istream& input,
                                    const std::string_view description) {
  std::string line;
  if (!std::getline(input, line)) {
    throw std::runtime_error("checkpoint is truncated before " +
                             std::string(description));
  }
  return line;
}

[[nodiscard]] std::string_view value_after_key(const std::string_view line,
                                               const std::string_view key) {
  const std::string prefix = std::string(key) + " = ";
  if (!line.starts_with(prefix)) {
    throw std::runtime_error("checkpoint expected key " + std::string(key));
  }
  const auto value = line.substr(prefix.size());
  if (value.empty()) {
    throw std::runtime_error("checkpoint value for " + std::string(key) + " is empty");
  }
  return value;
}

template <typename Value>
[[nodiscard]] Value parse_number(const std::string_view text,
                                 const std::string_view description) {
  Value value{};
  const auto result = std::from_chars(text.data(), text.data() + text.size(), value);
  if (result.ec != std::errc{} || result.ptr != text.data() + text.size()) {
    throw std::runtime_error("invalid checkpoint " + std::string(description));
  }
  return value;
}

}  // namespace

void write_checkpoint(std::ostream& output, const Checkpoint& checkpoint) {
  require_non_negative(checkpoint.time_s, "checkpoint time_s");
  if (checkpoint.state.empty()) {
    throw std::invalid_argument("checkpoint state must not be empty");
  }
  if (checkpoint.config_fingerprint.empty()) {
    throw std::invalid_argument("checkpoint config fingerprint must not be empty");
  }
  if (checkpoint.layout_id.empty() ||
      checkpoint.layout_id.find_first_of("\r\n") != std::string::npos) {
    throw std::invalid_argument("checkpoint layout ID must be a nonempty line");
  }
  for (const Real value : checkpoint.state) {
    require_finite(value, "checkpoint state value");
  }

  output.imbue(std::locale::classic());
  const bool legacy = checkpoint.layout_id == "legacy_scalar_v1";
  output << "MPS_CHECKPOINT\n"
         << "version = " << (legacy ? 1 : 2) << '\n'
         << "config_fingerprint = " << checkpoint.config_fingerprint << '\n';
  if (!legacy) {
    output << "layout_id = " << checkpoint.layout_id << '\n';
  }
  output << std::setprecision(std::numeric_limits<Real>::max_digits10)
         << "time_s = " << checkpoint.time_s << '\n'
         << "step = " << checkpoint.step << '\n'
         << "state_size = " << checkpoint.state.size() << '\n'
         << "state.begin\n";
  for (const Real value : checkpoint.state) {
    output << value << '\n';
  }
  output << "state.end\n";
  if (!output) {
    throw std::runtime_error("failed while writing checkpoint");
  }
}

Checkpoint read_checkpoint(std::istream& input,
                           const std::string_view expected_config_fingerprint,
                           const std::string_view expected_layout_id,
                           const std::size_t expected_state_size) {
  input.imbue(std::locale::classic());
  if (read_line(input, "magic") != "MPS_CHECKPOINT") {
    throw std::runtime_error("invalid checkpoint magic");
  }

  const auto version = parse_number<unsigned int>(
      value_after_key(read_line(input, "version"), "version"), "version");
  if (version != 1 && version != 2) {
    throw std::runtime_error("unsupported checkpoint version " +
                             std::to_string(version));
  }

  const std::string fingerprint(
      value_after_key(read_line(input, "config fingerprint"), "config_fingerprint"));
  if (fingerprint != expected_config_fingerprint) {
    throw std::runtime_error("checkpoint configuration fingerprint mismatch");
  }

  const std::string layout_id =
      version == 1
          ? "legacy_scalar_v1"
          : std::string(value_after_key(read_line(input, "layout ID"), "layout_id"));
  if (!expected_layout_id.empty() && layout_id != expected_layout_id) {
    throw std::runtime_error("checkpoint state layout mismatch");
  }

  const auto time_s =
      parse_number<Real>(value_after_key(read_line(input, "time"), "time_s"), "time_s");
  const auto step = parse_number<std::uint64_t>(
      value_after_key(read_line(input, "step"), "step"), "step");
  const auto state_size = parse_number<std::size_t>(
      value_after_key(read_line(input, "state size"), "state_size"), "state_size");
  require_non_negative(time_s, "checkpoint time_s");
  if (state_size == 0 || state_size > kMaximumStateSize) {
    throw std::runtime_error("checkpoint state size is outside the supported range");
  }
  if (expected_state_size != 0 && state_size != expected_state_size) {
    throw std::runtime_error("checkpoint state size does not match expected layout");
  }
  if (read_line(input, "state begin marker") != "state.begin") {
    throw std::runtime_error("checkpoint expected state.begin");
  }

  std::vector<Real> state;
  state.reserve(state_size);
  for (std::size_t index = 0; index < state_size; ++index) {
    const auto value =
        parse_number<Real>(read_line(input, "state value"), "state value");
    require_finite(value, "checkpoint state value");
    state.push_back(value);
  }
  if (read_line(input, "state end marker") != "state.end") {
    throw std::runtime_error("checkpoint expected state.end");
  }

  std::string trailing;
  if (std::getline(input, trailing)) {
    throw std::runtime_error("checkpoint contains trailing data");
  }
  if (input.bad()) {
    throw std::runtime_error("failed while reading checkpoint");
  }

  return Checkpoint{.time_s = time_s,
                    .step = step,
                    .state = std::move(state),
                    .config_fingerprint = fingerprint,
                    .layout_id = layout_id};
}

void write_checkpoint_file(const std::filesystem::path& path,
                           const Checkpoint& checkpoint) {
  auto temporary_path = path;
  temporary_path += ".tmp";
  try {
    std::ofstream output(temporary_path, std::ios::trunc);
    if (!output) {
      throw std::runtime_error("unable to open temporary checkpoint: " +
                               temporary_path.string());
    }
    write_checkpoint(output, checkpoint);
    output.close();
    if (!output) {
      throw std::runtime_error("unable to close temporary checkpoint: " +
                               temporary_path.string());
    }
    std::filesystem::rename(temporary_path, path);
  } catch (...) {
    std::error_code ignored;
    std::filesystem::remove(temporary_path, ignored);
    throw;
  }
}

Checkpoint read_checkpoint_file(const std::filesystem::path& path,
                                const std::string_view expected_config_fingerprint,
                                const std::string_view expected_layout_id,
                                const std::size_t expected_state_size) {
  std::ifstream input(path);
  if (!input) {
    throw std::runtime_error("unable to open checkpoint file: " + path.string());
  }
  return read_checkpoint(input, expected_config_fingerprint, expected_layout_id,
                         expected_state_size);
}

}  // namespace mps
