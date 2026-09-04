#pragma once

#include <optional>
#include <vector>

#include "myplanetsim/config/experiment_config.hpp"
#include "myplanetsim/dynamics/vertical_column_state.hpp"
#include "myplanetsim/vertical/hybrid_pressure_coordinate.hpp"
#include "myplanetsim/vertical/vertical_transport.hpp"

namespace mps {

struct VerticalColumnSample {
  Real time_s;
  std::uint64_t step;
  VerticalColumnState state;
  VerticalMassFlux mass_flux;
  Real maximum_cfl;
};

struct VerticalColumnResult {
  VerticalColumnState state;
  bool reached_end_time;
  Real maximum_cfl;
  std::vector<VerticalColumnSample> samples;
};

[[nodiscard]] AtmosphericHybridCoordinate make_vertical_coordinate(
    const ExperimentConfig& config);
[[nodiscard]] VerticalColumnState make_vertical_column_initial_state(
    const ExperimentConfig& config, const AtmosphericHybridCoordinate& coordinate);
[[nodiscard]] VerticalColumnResult run_vertical_column(
    const ExperimentConfig& config,
    std::optional<VerticalColumnState> initial_state = std::nullopt,
    std::optional<std::uint64_t> stop_after_step = std::nullopt);

}  // namespace mps
