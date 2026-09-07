#include "myplanetsim/io/run_metadata.hpp"

#include <cstdint>
#include <iomanip>
#include <locale>
#include <sstream>
#include <stdexcept>
#include <string>

#include "myplanetsim/build_info.hpp"
#include "myplanetsim/version.hpp"

namespace mps {
namespace {

[[nodiscard]] std::string hexadecimal(const std::uint64_t value) {
  std::ostringstream output;
  output << std::hex << std::setfill('0') << std::setw(16) << value;
  return output.str();
}

}  // namespace

std::string canonical_config_text(const ExperimentConfig& config) {
  std::ostringstream output;
  write_experiment_config(output, config);
  return output.str();
}

std::string config_fingerprint(const ExperimentConfig& config) {
  constexpr std::uint64_t kFnvOffsetBasis = 14695981039346656037ULL;
  constexpr std::uint64_t kFnvPrime = 1099511628211ULL;
  std::uint64_t hash = kFnvOffsetBasis;
  for (const unsigned char byte : canonical_config_text(config)) {
    hash ^= byte;
    hash *= kFnvPrime;
  }
  return hexadecimal(hash);
}

RunMetadata make_run_metadata(const ExperimentConfig& config) {
  config.validate();
  return RunMetadata{
      .model_version = std::string(version()),
      .git_revision = std::string(build_info::git_revision),
      .git_dirty = build_info::git_dirty,
      .compiler_id = std::string(build_info::compiler_id),
      .compiler_version = std::string(build_info::compiler_version),
      .build_type = std::string(build_info::build_type),
      .random_seed = config.run.random_seed,
      .config_fingerprint = config_fingerprint(config),
  };
}

void write_run_metadata(std::ostream& output, const RunMetadata& metadata,
                        const ExperimentConfig& config) {
  output.imbue(std::locale::classic());
  output << "metadata.format_version = 1\n"
         << "metadata.model_version = " << metadata.model_version << '\n'
         << "metadata.git_revision = " << metadata.git_revision << '\n'
         << "metadata.git_dirty = " << std::boolalpha << metadata.git_dirty << '\n'
         << "metadata.compiler_id = " << metadata.compiler_id << '\n'
         << "metadata.compiler_version = " << metadata.compiler_version << '\n'
         << "metadata.build_type = " << metadata.build_type << '\n'
         << "metadata.random_seed = " << metadata.random_seed << '\n'
         << "metadata.config_fingerprint = " << metadata.config_fingerprint << '\n';
  if (!config.tracers.empty()) {
    output << "metadata.tracer_count = " << config.tracers.size() << '\n';
    for (std::size_t index = 0; index < config.tracers.size(); ++index)
      output << "metadata.tracer." << index << ".name = " << config.tracers[index].name
             << '\n'
             << "metadata.tracer." << index
             << ".role = " << tracer_role_name(config.tracers[index].role) << '\n';
  }
  output << "configuration.begin\n"
         << canonical_config_text(config) << "configuration.end\n";
  if (!output) {
    throw std::runtime_error("failed while writing run metadata");
  }
}

}  // namespace mps
