#pragma once

#include <cstdint>
#include <filesystem>
#include <iosfwd>
#include <string>

#include "myplanetsim/config/experiment_config.hpp"

namespace mps {

inline constexpr std::uint64_t kControlProtocolVersion = 1;
inline constexpr std::uint64_t kControlMaxEditCount = 64;

struct ControlRequestV1 {
  std::uint64_t format_version = kControlProtocolVersion;
  std::string run_id;
  Real end_time_s = 0.0;
  Real maximum_time_step_s = 0.0;
  std::uint64_t frame_interval_steps = 0;
  std::filesystem::path frame_directory;
  std::uint64_t initial_edits_count = 0;
};

[[nodiscard]] ControlRequestV1 parse_control_request(std::istream& input);
[[nodiscard]] ControlRequestV1 load_control_request(const std::filesystem::path& path);
void write_control_request(std::ostream& output, const ControlRequestV1& request);

void validate_control_request(const ControlRequestV1& request);
[[nodiscard]] ExperimentConfig apply_control_request(const ExperimentConfig& config,
                                                     const ControlRequestV1& request);

}  // namespace mps
