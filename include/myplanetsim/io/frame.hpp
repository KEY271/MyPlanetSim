#pragma once

#include <cstdint>
#include <filesystem>
#include <string>
#include <string_view>

#include "myplanetsim/dynamics/dry_hydrostatic_state.hpp"
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

struct FrameV2Metadata {
  std::uint32_t schema_version = 2;
  Index cells_per_panel = 0;
  Index levels = 0;
  Real time_s = 0;
  std::uint64_t step = 0;
  std::string config_fingerprint;
};
struct FrameV2 {
  FrameV2Metadata metadata;
  std::vector<Real> surface_pressure_pa;
  DryHydrostaticDerived derived;
  std::uint64_t byte_length = 0;
};

void write_frame_file(const std::filesystem::path& path,
                      const FrameV1Metadata& metadata, const ShallowWaterState& state);
[[nodiscard]] FrameV1 read_frame_file(
    const std::filesystem::path& path,
    std::string_view expected_config_fingerprint = {});
void write_frame_v2_file(const std::filesystem::path&, const FrameV2Metadata&,
                         std::span<const Real> surface_pressure_pa,
                         const DryHydrostaticDerived&);
[[nodiscard]] FrameV2 read_frame_v2_file(
    const std::filesystem::path&, std::string_view expected_config_fingerprint = {});

}  // namespace mps
