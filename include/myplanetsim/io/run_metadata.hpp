#pragma once

#include <iosfwd>
#include <string>

#include "myplanetsim/config/experiment_config.hpp"
#include "myplanetsim/core/types.hpp"

namespace mps {

struct RunMetadata {
  std::string model_version;
  std::string git_revision;
  bool git_dirty;
  std::string compiler_id;
  std::string compiler_version;
  std::string build_type;
  Seed random_seed;
  std::string config_fingerprint;
};

[[nodiscard]] std::string canonical_config_text(const ExperimentConfig& config);
[[nodiscard]] std::string config_fingerprint(const ExperimentConfig& config);
[[nodiscard]] RunMetadata make_run_metadata(const ExperimentConfig& config);
void write_run_metadata(std::ostream& output, const RunMetadata& metadata,
                        const ExperimentConfig& config);
}  // namespace mps
