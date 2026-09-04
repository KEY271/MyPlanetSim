#pragma once

#include <cstdint>
#include <filesystem>
#include <iosfwd>
#include <string>
#include <vector>

#include "myplanetsim/config/experiment_config.hpp"
#include "myplanetsim/geometry/vec3.hpp"

namespace mps {

inline constexpr std::uint64_t kControlProtocolVersion = 1;
inline constexpr std::uint64_t kControlMaxEditCount = 64;
inline constexpr std::uint64_t kControlMaxCellsPerPanel = 96;
inline constexpr std::uint64_t kControlMaxPublishedFrames = 512;

enum class InitialConditionEditKind { kGaussianDepth };
enum class MassPolicy { kPreserveGlobal, kAllowChange };

struct InitialConditionEditV1 {
  InitialConditionEditKind kind = InitialConditionEditKind::kGaussianDepth;
  Vec3 center_unit{0.0, 0.0, 1.0};
  Real amplitude_m = 0.0;
  Real sigma_rad = 0.0;
  MassPolicy mass_policy = MassPolicy::kPreserveGlobal;
};

struct ControlRequestV1 {
  std::uint64_t format_version = kControlProtocolVersion;
  std::string run_id;
  std::uint64_t cells_per_panel = 0;
  Real end_time_s = 0.0;
  Real maximum_time_step_s = 0.0;
  std::uint64_t frame_interval_steps = 0;
  std::filesystem::path frame_directory;
  std::vector<InitialConditionEditV1> initial_edits;
};

[[nodiscard]] ControlRequestV1 parse_control_request(std::istream& input);
[[nodiscard]] ControlRequestV1 load_control_request(const std::filesystem::path& path);
void write_control_request(std::ostream& output, const ControlRequestV1& request);

void validate_control_request(const ControlRequestV1& request);
[[nodiscard]] ExperimentConfig apply_control_request(const ExperimentConfig& config,
                                                     const ControlRequestV1& request);

[[nodiscard]] std::string_view initial_condition_edit_kind_name(
    InitialConditionEditKind kind) noexcept;
[[nodiscard]] std::string_view mass_policy_name(MassPolicy policy) noexcept;

}  // namespace mps
