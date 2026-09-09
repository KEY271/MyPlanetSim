#include "myplanetsim/diagnostics/period_mean.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "myplanetsim/core/validation.hpp"

namespace mps::diagnostics {

PeriodMeanAccumulator::PeriodMeanAccumulator(const std::size_t value_count,
                                             const Real start_time_s,
                                             const Real period_s)
    : value_count_(value_count),
      start_time_s_(start_time_s),
      period_s_(period_s),
      last_step_end_s_(start_time_s),
      integral_(value_count) {
  if (value_count == 0) throw std::invalid_argument("period mean requires values");
  require_finite(start_time_s, "period mean start time");
  require_positive(period_s, "period mean period");
}

PeriodMeanAccumulator::PeriodMeanAccumulator(PeriodMeanAccumulatorState state)
    : value_count_(state.value_count),
      start_time_s_(state.start_time_s),
      period_s_(state.period_s),
      last_step_end_s_(state.last_step_end_s),
      has_observation_(state.has_observation),
      open_index_(state.open_index),
      open_actual_start_s_(state.open_actual_start_s),
      open_actual_end_s_(state.open_actual_end_s),
      open_weight_s_(state.open_weight_s),
      integral_(std::move(state.integral)) {
  if (value_count_ == 0 || integral_.size() != value_count_)
    throw std::invalid_argument("invalid period mean restart shape");
  require_finite(start_time_s_, "period mean restart start time");
  require_positive(period_s_, "period mean restart period");
  require_non_negative(open_weight_s_, "period mean restart weight");
  if (has_observation_) require_finite(last_step_end_s_, "period mean last step end");
  if ((open_weight_s_ == 0.0) != (open_actual_end_s_ == 0.0))
    throw std::invalid_argument("invalid period mean restart bounds");
  for (const auto value : integral_) require_finite(value, "period mean integral");
}

PeriodMeanWindow PeriodMeanAccumulator::make_window(const bool complete) const {
  PeriodMeanWindow result{
      .index = open_index_,
      .scheduled_start_s = start_time_s_ + static_cast<Real>(open_index_) * period_s_,
      .scheduled_end_s = start_time_s_ + static_cast<Real>(open_index_ + 1) * period_s_,
      .actual_start_s = open_actual_start_s_,
      .actual_end_s = open_actual_end_s_,
      .weight_s = open_weight_s_,
      .complete = complete,
      .values = integral_,
  };
  for (auto& value : result.values) value /= open_weight_s_;
  return result;
}

void PeriodMeanAccumulator::close_open_window() {
  if (open_weight_s_ > 0.0) completed_.push_back(make_window(true));
  ++open_index_;
  open_weight_s_ = 0.0;
  open_actual_start_s_ = 0.0;
  open_actual_end_s_ = 0.0;
  std::ranges::fill(integral_, 0.0);
}

void PeriodMeanAccumulator::observe(const Real step_start_s, const Real step_end_s,
                                    const std::span<const Real> values) {
  require_finite(step_start_s, "period mean step start");
  require_finite(step_end_s, "period mean step end");
  if (values.size() != value_count_)
    throw std::invalid_argument("period mean value count mismatch");
  if (step_end_s < step_start_s)
    throw std::invalid_argument("period mean step end precedes start");
  if (has_observation_ && step_start_s < last_step_end_s_)
    throw std::invalid_argument("period mean observations overlap or go backwards");
  for (const Real value : values)
    if (!std::isfinite(value))
      throw std::invalid_argument("period mean values must be finite");
  if (step_end_s == step_start_s) return;
  has_observation_ = true;
  last_step_end_s_ = step_end_s;

  Real cursor = std::max(step_start_s, start_time_s_);
  if (step_end_s <= cursor) return;
  open_index_ =
      static_cast<std::uint64_t>(std::floor((cursor - start_time_s_) / period_s_));
  while (cursor < step_end_s) {
    const Real window_start =
        start_time_s_ + static_cast<Real>(open_index_) * period_s_;
    const Real window_end = window_start + period_s_;
    // Exact boundary values belong to the next half-open window.
    if (cursor >= window_end) {
      close_open_window();
      continue;
    }
    const Real overlap_end = std::min(step_end_s, window_end);
    const Real weight = overlap_end - cursor;
    if (weight > 0.0) {
      if (open_weight_s_ == 0.0) open_actual_start_s_ = cursor;
      open_actual_end_s_ = overlap_end;
      open_weight_s_ += weight;
      for (std::size_t i = 0; i < value_count_; ++i) integral_[i] += weight * values[i];
    }
    cursor = overlap_end;
    if (cursor >= window_end) close_open_window();
  }
}

std::vector<PeriodMeanWindow> PeriodMeanAccumulator::take_completed() {
  auto result = std::move(completed_);
  completed_.clear();
  return result;
}

std::optional<PeriodMeanWindow> PeriodMeanAccumulator::partial() const {
  if (open_weight_s_ == 0.0) return std::nullopt;
  return make_window(false);
}

PeriodMeanAccumulatorState PeriodMeanAccumulator::state() const {
  return {.value_count = value_count_,
          .start_time_s = start_time_s_,
          .period_s = period_s_,
          .last_step_end_s = last_step_end_s_,
          .has_observation = has_observation_,
          .open_index = open_index_,
          .open_actual_start_s = open_actual_start_s_,
          .open_actual_end_s = open_actual_end_s_,
          .open_weight_s = open_weight_s_,
          .integral = integral_};
}

}  // namespace mps::diagnostics
