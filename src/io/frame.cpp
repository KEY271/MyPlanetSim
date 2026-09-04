#include "myplanetsim/io/frame.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstring>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <vector>

#include "myplanetsim/core/validation.hpp"

namespace mps {
namespace {

constexpr std::array<unsigned char, 8> kMagic{'M', 'P', 'S', 'F', 'R', 'A', 'M', '1'};
constexpr std::uint32_t kVersion = 1;
constexpr std::uint64_t kMaximumPayloadBytes = 1024ULL * 1024ULL * 1024ULL;

void append_u32(std::vector<unsigned char>& bytes, const std::uint32_t value) {
  for (unsigned int shift = 0; shift < 32; shift += 8) {
    bytes.push_back(static_cast<unsigned char>(value >> shift));
  }
}

void append_u64(std::vector<unsigned char>& bytes, const std::uint64_t value) {
  for (unsigned int shift = 0; shift < 64; shift += 8) {
    bytes.push_back(static_cast<unsigned char>(value >> shift));
  }
}

[[nodiscard]] std::uint32_t read_u32(const std::vector<unsigned char>& bytes,
                                     std::size_t& offset) {
  if (offset + 4 > bytes.size()) {
    throw std::runtime_error("frame is truncated");
  }
  std::uint32_t value = 0;
  for (unsigned int shift = 0; shift < 32; shift += 8) {
    value |= static_cast<std::uint32_t>(bytes[offset++]) << shift;
  }
  return value;
}

[[nodiscard]] std::uint64_t read_u64(const std::vector<unsigned char>& bytes,
                                     std::size_t& offset) {
  if (offset + 8 > bytes.size()) {
    throw std::runtime_error("frame is truncated");
  }
  std::uint64_t value = 0;
  for (unsigned int shift = 0; shift < 64; shift += 8) {
    value |= static_cast<std::uint64_t>(bytes[offset++]) << shift;
  }
  return value;
}

void append_real(std::vector<unsigned char>& bytes, const Real value) {
  append_u64(bytes, std::bit_cast<std::uint64_t>(value));
}

[[nodiscard]] Real read_real(const std::vector<unsigned char>& bytes,
                             std::size_t& offset) {
  return std::bit_cast<Real>(read_u64(bytes, offset));
}

[[nodiscard]] std::vector<unsigned char> read_bytes(const std::filesystem::path& path) {
  std::ifstream input(path, std::ios::binary | std::ios::ate);
  if (!input) {
    throw std::runtime_error("unable to open frame: " + path.string());
  }
  const auto size = input.tellg();
  if (size < 0 || static_cast<std::uint64_t>(size) > kMaximumPayloadBytes) {
    throw std::runtime_error("frame is outside the supported size");
  }
  std::vector<unsigned char> bytes(static_cast<std::size_t>(size));
  input.seekg(0);
  if (!bytes.empty()) {
    input.read(reinterpret_cast<char*>(bytes.data()),
               static_cast<std::streamsize>(bytes.size()));
  }
  if (!input) {
    throw std::runtime_error("failed while reading frame");
  }
  return bytes;
}

}  // namespace

void write_frame_file(const std::filesystem::path& path,
                      const FrameV1Metadata& metadata, const ShallowWaterState& state) {
  if (metadata.schema_version != kVersion || metadata.cells_per_panel <= 0 ||
      metadata.cells_per_panel > 1024 || metadata.config_fingerprint.empty() ||
      metadata.config_fingerprint.size() > 128) {
    throw std::invalid_argument("invalid FrameV1 metadata");
  }
  require_non_negative(metadata.time_s, "frame time_s");
  if (state.depth.empty() || state.momentum.size() != state.depth.size()) {
    throw std::invalid_argument("frame state shape is invalid");
  }
  const auto n = static_cast<std::uint64_t>(metadata.cells_per_panel);
  const auto expected_cells = 6 * n * n;
  if (state.depth.size() != expected_cells) {
    throw std::invalid_argument("frame state size does not match cells_per_panel");
  }
  std::vector<unsigned char> bytes;
  bytes.reserve(64 + metadata.config_fingerprint.size() + 32 * state.depth.size());
  bytes.insert(bytes.end(), kMagic.begin(), kMagic.end());
  append_u32(bytes, kVersion);
  append_u32(bytes, 0);
  append_u64(bytes, n);
  append_real(bytes, metadata.time_s);
  append_u64(bytes, metadata.step);
  append_u64(bytes, state.depth.size());
  append_u32(bytes, static_cast<std::uint32_t>(metadata.config_fingerprint.size()));
  bytes.insert(bytes.end(), metadata.config_fingerprint.begin(),
               metadata.config_fingerprint.end());
  const auto flat = flatten_shallow_water_state(state);
  for (const Real value : flat) {
    require_finite(value, "frame value");
    append_real(bytes, value);
  }

  if (!path.parent_path().empty()) {
    std::filesystem::create_directories(path.parent_path());
  }
  auto temporary_path = path;
  temporary_path += ".tmp";
  std::ofstream output(temporary_path, std::ios::binary | std::ios::trunc);
  if (!output) {
    throw std::runtime_error("unable to open temporary frame: " +
                             temporary_path.string());
  }
  output.write(reinterpret_cast<const char*>(bytes.data()),
               static_cast<std::streamsize>(bytes.size()));
  output.close();
  if (!output) {
    throw std::runtime_error("failed while writing frame");
  }
  std::filesystem::rename(temporary_path, path);
}

FrameV1 read_frame_file(const std::filesystem::path& path,
                        const std::string_view expected_config_fingerprint) {
  const auto bytes = read_bytes(path);
  std::size_t offset = 0;
  if (bytes.size() < kMagic.size() ||
      !std::equal(kMagic.begin(), kMagic.end(), bytes.begin())) {
    throw std::runtime_error("invalid FrameV1 magic");
  }
  offset = kMagic.size();
  if (read_u32(bytes, offset) != kVersion || read_u32(bytes, offset) != 0) {
    throw std::runtime_error("unsupported or invalid FrameV1 header");
  }
  const auto cells_per_panel = read_u64(bytes, offset);
  if (cells_per_panel == 0 || cells_per_panel > 1024) {
    throw std::runtime_error("FrameV1 cells_per_panel is outside the supported range");
  }
  const Real time_s = read_real(bytes, offset);
  const auto step = read_u64(bytes, offset);
  const auto cell_count = read_u64(bytes, offset);
  const auto fingerprint_size = read_u32(bytes, offset);
  if (fingerprint_size == 0 || fingerprint_size > 128 ||
      offset + fingerprint_size > bytes.size()) {
    throw std::runtime_error("invalid FrameV1 fingerprint");
  }
  const std::string fingerprint(reinterpret_cast<const char*>(bytes.data() + offset),
                                fingerprint_size);
  offset += fingerprint_size;
  if (!expected_config_fingerprint.empty() &&
      fingerprint != expected_config_fingerprint) {
    throw std::runtime_error("frame configuration fingerprint mismatch");
  }
  const auto expected_cells = 6 * cells_per_panel * cells_per_panel;
  if (cell_count != expected_cells ||
      cell_count > std::numeric_limits<std::size_t>::max() / 4) {
    throw std::runtime_error("FrameV1 cell count does not match shape");
  }
  const auto value_count = static_cast<std::size_t>(4 * cell_count);
  if (offset + value_count * sizeof(Real) != bytes.size()) {
    throw std::runtime_error("FrameV1 payload has the wrong size");
  }
  std::vector<Real> flat;
  flat.reserve(value_count);
  for (std::size_t index = 0; index < value_count; ++index) {
    const Real value = read_real(bytes, offset);
    require_finite(value, "frame value");
    flat.push_back(value);
  }
  return {.metadata = {.schema_version = kVersion,
                       .cells_per_panel = static_cast<Index>(cells_per_panel),
                       .time_s = time_s,
                       .step = step,
                       .config_fingerprint = fingerprint},
          .state = unflatten_shallow_water_state(time_s, step, flat,
                                                 static_cast<std::size_t>(cell_count)),
          .byte_length = static_cast<std::uint64_t>(bytes.size())};
}

}  // namespace mps
