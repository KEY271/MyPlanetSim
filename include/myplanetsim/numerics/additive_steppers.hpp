#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

#include "myplanetsim/numerics/explicit_steppers.hpp"

namespace mps {

struct ThreeStageAdditiveButcherTable {
  std::array<Real, 3> nodes{};
  std::array<std::array<Real, 3>, 3> explicit_matrix{};
  std::array<std::array<Real, 3>, 3> implicit_matrix{};
  std::array<Real, 3> explicit_weights{};
  std::array<Real, 3> implicit_weights{};
  std::size_t order = 0;
};

// SUNDIALS ARKODE v7.7.0 ARK2 ERK/DIRK pair. Keeping the complete table here
// makes the stage times and split weights reviewable without a runtime dependency.
[[nodiscard]] inline ThreeStageAdditiveButcherTable ark2_imex_table() {
  const Real root_two = std::sqrt(2.0);
  const Real gamma = 1.0 - 1.0 / root_two;
  const Real diagonal_weight = 1.0 / (2.0 * root_two);
  const Real explicit_third_stage_second = (3.0 + 2.0 * root_two) / 6.0;
  return {.nodes = {0.0, 2.0 - root_two, 1.0},
          .explicit_matrix = {{{0.0, 0.0, 0.0},
                               {2.0 - root_two, 0.0, 0.0},
                               {1.0 - explicit_third_stage_second,
                                explicit_third_stage_second, 0.0}}},
          .implicit_matrix = {{{0.0, 0.0, 0.0},
                               {gamma, gamma, 0.0},
                               {diagonal_weight, diagonal_weight, gamma}}},
          .explicit_weights = {diagonal_weight, diagonal_weight, gamma},
          .implicit_weights = {diagonal_weight, diagonal_weight, gamma},
          .order = 2};
}

// Fixed-step comparator for y' = E(t,y) + I(t,y), where I is linear. The
// caller-provided solve computes (identity - implicit_time_s * I) stage = rhs.
// ARK2 has three explicit/implicit evaluations and two non-trivial stage solves.
class Ark2Imex {
 public:
  explicit Ark2Imex(const std::size_t state_size) : right_hand_side_(state_size) {
    if (state_size == 0)
      throw std::invalid_argument("stepper state size must be positive");
    for (auto& value : stages_) value.resize(state_size);
    for (auto& value : explicit_tendencies_) value.resize(state_size);
    for (auto& value : implicit_tendencies_) value.resize(state_size);
  }

  [[nodiscard]] std::size_t state_size() const noexcept {
    return right_hand_side_.size();
  }

  template <typename ExplicitRhs, typename ImplicitRhs, typename ImplicitSolve>
  void step(const Real time_s, const Real time_step_s, const std::span<Real> state,
            ExplicitRhs&& explicit_rhs, ImplicitRhs&& implicit_rhs,
            ImplicitSolve&& implicit_solve) {
    detail::validate_step_arguments(time_s, time_step_s, state, state_size());
    const auto table = ark2_imex_table();
    std::copy(state.begin(), state.end(), stages_[0].begin());
    evaluate_stage(0, time_s, explicit_rhs, implicit_rhs);

    for (std::size_t stage = 1; stage < stages_.size(); ++stage) {
      std::copy(stages_[0].begin(), stages_[0].end(), right_hand_side_.begin());
      for (std::size_t previous = 0; previous < stage; ++previous) {
        const Real explicit_coefficient = table.explicit_matrix[stage][previous];
        const Real implicit_coefficient = table.implicit_matrix[stage][previous];
        for (std::size_t index = 0; index < state_size(); ++index) {
          right_hand_side_[index] +=
              time_step_s *
              (explicit_coefficient * explicit_tendencies_[previous][index] +
               implicit_coefficient * implicit_tendencies_[previous][index]);
        }
      }
      const Real stage_time = time_s + table.nodes[stage] * time_step_s;
      implicit_solve(stage_time, time_step_s * table.implicit_matrix[stage][stage],
                     std::span<const Real>(right_hand_side_),
                     std::span<Real>(stages_[stage]));
      detail::validate_tendency(stages_[stage]);
      evaluate_stage(stage, stage_time, explicit_rhs, implicit_rhs);
    }

    for (std::size_t index = 0; index < state_size(); ++index) {
      Real next = stages_[0][index];
      for (std::size_t stage = 0; stage < stages_.size(); ++stage) {
        next += time_step_s *
                (table.explicit_weights[stage] * explicit_tendencies_[stage][index] +
                 table.implicit_weights[stage] * implicit_tendencies_[stage][index]);
      }
      state[index] = next;
    }
    detail::validate_tendency(state);
  }

 private:
  template <typename ExplicitRhs, typename ImplicitRhs>
  void evaluate_stage(const std::size_t stage, const Real time_s,
                      ExplicitRhs& explicit_rhs, ImplicitRhs& implicit_rhs) {
    explicit_rhs(time_s, std::span<const Real>(stages_[stage]),
                 std::span<Real>(explicit_tendencies_[stage]));
    detail::validate_tendency(explicit_tendencies_[stage]);
    implicit_rhs(time_s, std::span<const Real>(stages_[stage]),
                 std::span<Real>(implicit_tendencies_[stage]));
    detail::validate_tendency(implicit_tendencies_[stage]);
  }

  std::array<std::vector<Real>, 3> stages_;
  std::array<std::vector<Real>, 3> explicit_tendencies_;
  std::array<std::vector<Real>, 3> implicit_tendencies_;
  std::vector<Real> right_hand_side_;
};

}  // namespace mps
