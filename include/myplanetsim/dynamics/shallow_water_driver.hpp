#pragma once

#include <cstdint>
#include <optional>
#include <span>
#include <vector>

#include "myplanetsim/config/experiment_config.hpp"
#include "myplanetsim/control/control_request.hpp"
#include "myplanetsim/diagnostics/shallow_water_diagnostics.hpp"

namespace mps {

struct ShallowWaterSample {
  Real time_s;
  std::uint64_t step;
  diagnostics::ShallowWaterInvariants invariants;
};

struct ShallowWaterResult {
  ShallowWaterState state;
  diagnostics::ShallowWaterInvariants initial_diagnostics;
  diagnostics::ShallowWaterInvariants final_diagnostics;
  bool reached_end_time;
  Real maximum_cfl;
  // Invariants sampled every diagnostics.interval_steps, ends included.
  std::vector<ShallowWaterSample> samples;
};

struct ShallowWaterRunHooks {
  using FrameObserver = void (*)(const ShallowWaterState&, void* context);
  using CancellationPredicate = bool (*)(void* context);

  FrameObserver on_frame = nullptr;
  void* observer_context = nullptr;
  CancellationPredicate is_cancelled = nullptr;
  void* cancellation_context = nullptr;
};

[[nodiscard]] Real stable_shallow_water_time_step(const CubedSphereGrid& grid,
                                                  const ShallowWaterState& state,
                                                  Real gravity_m_s2, Real cfl,
                                                  Real maximum_time_step_s);
[[nodiscard]] Real shallow_water_cfl_number(const CubedSphereGrid& grid,
                                            const ShallowWaterState& state,
                                            Real gravity_m_s2, Real time_step_s);
[[nodiscard]] ShallowWaterResult run_shallow_water(
    const ExperimentConfig& config,
    std::optional<ShallowWaterState> initial_state = std::nullopt,
    std::optional<std::uint64_t> stop_after_step = std::nullopt,
    std::span<const InitialConditionEditV1> initial_edits = {},
    ShallowWaterRunHooks hooks = {});

}  // namespace mps
