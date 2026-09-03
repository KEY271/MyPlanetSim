#pragma once

#include <cstdint>
#include <filesystem>
#include <string>
#include <string_view>

#include "myplanetsim/dynamics/shallow_water_state.hpp"

namespace mps {

struct FrameV1Metadata {
  std::uint32_t schema_version = 1;
  Index cells_per_panel = 0;
  Real time_s = 0.0;
  std::uint64_t step = 0;
  std::string config_fingerprint;
};

struct FrameV1 {
  FrameV1Metadata metadata;
  ShallowWaterState state;
  std::uint64_t byte_length = 0;
};

void write_frame_file(const std::filesystem::path& path,
                      const FrameV1Metadata& metadata,
                      const ShallowWaterState& state);
[[nodiscard]] FrameV1 read_frame_file(
    const std::filesystem::path& path,
    std::string_view expected_config_fingerprint = {});

}  // namespace mps
