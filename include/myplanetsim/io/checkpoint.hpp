#pragma once

#include <cstdint>
#include <filesystem>
#include <iosfwd>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#include "myplanetsim/core/types.hpp"

namespace mps {

struct Checkpoint {
  Real time_s;
  std::uint64_t step;
  std::vector<Real> state;
  std::string config_fingerprint;
  std::string layout_id = "legacy_scalar_v1";
};

void write_checkpoint(std::ostream& output, const Checkpoint& checkpoint);
[[nodiscard]] Checkpoint read_checkpoint(std::istream& input,
                                         std::string_view expected_config_fingerprint,
                                         std::string_view expected_layout_id = {},
                                         std::size_t expected_state_size = 0);
void write_checkpoint_file(const std::filesystem::path& path,
                           const Checkpoint& checkpoint);
[[nodiscard]] Checkpoint read_checkpoint_file(
    const std::filesystem::path& path, std::string_view expected_config_fingerprint,
    std::string_view expected_layout_id = {}, std::size_t expected_state_size = 0);

}  // namespace mps
