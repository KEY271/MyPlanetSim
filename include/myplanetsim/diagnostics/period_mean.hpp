#pragma once

#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <vector>

#include "myplanetsim/core/types.hpp"

namespace mps::diagnostics {

struct PeriodMeanWindow {
  std::uint64_t index = 0;
  Real scheduled_start_s = 0.0;
  Real scheduled_end_s = 0.0;
  Real actual_start_s = 0.0;
  Real actual_end_s = 0.0;
  Real weight_s = 0.0;
  bool complete = false;
  std::vector<Real> values;
};

struct PeriodMeanAccumulatorState {
  std::size_t value_count = 0;
  Real start_time_s = 0.0;
  Real period_s = 0.0;
  Real last_step_end_s = 0.0;
  bool has_observation = false;
  std::uint64_t open_index = 0;
  Real open_actual_start_s = 0.0;
  Real open_actual_end_s = 0.0;
  Real open_weight_s = 0.0;
  std::vector<Real> integral;
};

// Accumulates flattened fields using accepted-step end values as rectangle-rule
// representatives. One instance owns only the open window plus completed results.
class PeriodMeanAccumulator {
 public:
  PeriodMeanAccumulator(std::size_t value_count, Real start_time_s, Real period_s);
  explicit PeriodMeanAccumulator(PeriodMeanAccumulatorState state);

  void observe(Real step_start_s, Real step_end_s, std::span<const Real> values);
  [[nodiscard]] std::vector<PeriodMeanWindow> take_completed();
  [[nodiscard]] std::optional<PeriodMeanWindow> partial() const;
  [[nodiscard]] std::size_t value_count() const noexcept { return value_count_; }
  [[nodiscard]] Real start_time_s() const noexcept { return start_time_s_; }
  [[nodiscard]] Real period_s() const noexcept { return period_s_; }
  [[nodiscard]] PeriodMeanAccumulatorState state() const;

 private:
  std::size_t value_count_;
  Real start_time_s_;
  Real period_s_;
  Real last_step_end_s_;
  bool has_observation_ = false;
  std::uint64_t open_index_ = 0;
  Real open_actual_start_s_ = 0.0;
  Real open_actual_end_s_ = 0.0;
  Real open_weight_s_ = 0.0;
  std::vector<Real> integral_;
  std::vector<PeriodMeanWindow> completed_;

  [[nodiscard]] PeriodMeanWindow make_window(bool complete) const;
  void close_open_window();
};

}  // namespace mps::diagnostics
