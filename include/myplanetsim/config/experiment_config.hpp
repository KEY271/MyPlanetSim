#pragma once

#include <filesystem>
#include <iosfwd>
#include <string>

#include "myplanetsim/core/planet_parameters.hpp"
#include "myplanetsim/core/types.hpp"

namespace mps {

struct RunParameters {
  Real start_time_s;
  Real end_time_s;
  Real time_step_s;
  Seed random_seed;
};

struct OdeParameters {
  Real initial_value;
  Real decay_rate_s_1;
};

struct ExperimentConfig {
  PlanetParameters planet;
  RunParameters run;
  OdeParameters ode;
  std::string output_directory;

  void validate() const;
};

[[nodiscard]] ExperimentConfig parse_experiment_config(std::istream& input);
[[nodiscard]] ExperimentConfig load_experiment_config(
    const std::filesystem::path& path);
void write_experiment_config(std::ostream& output, const ExperimentConfig& config);

}  // namespace mps
