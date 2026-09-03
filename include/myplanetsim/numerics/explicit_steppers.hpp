#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <span>
#include <stdexcept>
#include <string_view>
#include <utility>
#include <vector>

#include "myplanetsim/core/types.hpp"
#include "myplanetsim/core/validation.hpp"

namespace mps {
namespace detail {

inline void validate_step_arguments(const Real time_s, const Real time_step_s,
                                    const std::span<const Real> state,
                                    const std::size_t expected_size) {
  require_finite(time_s, "time_s");
  require_positive(time_step_s, "time_step_s");
  if (state.size() != expected_size) {
    throw std::invalid_argument("state size does not match stepper size");
  }
  if (!std::ranges::all_of(state,
                           [](const Real value) { return std::isfinite(value); })) {
    throw std::invalid_argument("state contains a non-finite value");
  }
}

inline void validate_tendency(const std::span<const Real> tendency) {
  if (!std::ranges::all_of(tendency,
                           [](const Real value) { return std::isfinite(value); })) {
    throw std::runtime_error("RHS produced a non-finite tendency");
  }
}

}  // namespace detail

class ForwardEuler {
 public:
  explicit ForwardEuler(const std::size_t state_size) : tendency_(state_size) {
    if (state_size == 0) {
      throw std::invalid_argument("stepper state size must be positive");
    }
  }

  [[nodiscard]] std::size_t state_size() const noexcept { return tendency_.size(); }

  template <typename Rhs>
  void step(const Real time_s, const Real time_step_s, const std::span<Real> state,
            Rhs&& rhs) {
    detail::validate_step_arguments(time_s, time_step_s, state, state_size());
    std::forward<Rhs>(rhs)(time_s, std::span<const Real>(state),
                           std::span<Real>(tendency_));
    detail::validate_tendency(tendency_);
    for (std::size_t index = 0; index < state.size(); ++index) {
      state[index] += time_step_s * tendency_[index];
    }
  }

 private:
  std::vector<Real> tendency_;
};

class SspRk3 {
 public:
  explicit SspRk3(const std::size_t state_size)
      : initial_(state_size),
        stage_(state_size),
        next_stage_(state_size),
        tendency_(state_size) {
    if (state_size == 0) {
      throw std::invalid_argument("stepper state size must be positive");
    }
  }

  [[nodiscard]] std::size_t state_size() const noexcept { return tendency_.size(); }

  template <typename Rhs>
  void step(const Real time_s, const Real time_step_s, const std::span<Real> state,
            Rhs&& rhs) {
    detail::validate_step_arguments(time_s, time_step_s, state, state_size());
    std::copy(state.begin(), state.end(), initial_.begin());

    rhs(time_s, std::span<const Real>(initial_), std::span<Real>(tendency_));
    detail::validate_tendency(tendency_);
    for (std::size_t index = 0; index < state.size(); ++index) {
      stage_[index] = initial_[index] + time_step_s * tendency_[index];
    }

    rhs(time_s + time_step_s, std::span<const Real>(stage_),
        std::span<Real>(tendency_));
    detail::validate_tendency(tendency_);
    for (std::size_t index = 0; index < state.size(); ++index) {
      next_stage_[index] = 0.75 * initial_[index] +
                           0.25 * (stage_[index] + time_step_s * tendency_[index]);
    }

    rhs(time_s + 0.5 * time_step_s, std::span<const Real>(next_stage_),
        std::span<Real>(tendency_));
    detail::validate_tendency(tendency_);
    for (std::size_t index = 0; index < state.size(); ++index) {
      state[index] =
          (1.0 / 3.0) * initial_[index] +
          (2.0 / 3.0) * (next_stage_[index] + time_step_s * tendency_[index]);
    }
  }

 private:
  std::vector<Real> initial_;
  std::vector<Real> stage_;
  std::vector<Real> next_stage_;
  std::vector<Real> tendency_;
};

}  // namespace mps
