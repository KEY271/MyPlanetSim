#include <array>
#include <cmath>
#include <span>
#include <stdexcept>
#include <vector>

#include "myplanetsim/numerics/additive_steppers.hpp"
#include "support/test.hpp"

namespace {

struct SplitScalarProblem {
  mps::Real explicit_rate = 0.0;
  mps::Real implicit_rate = 0.0;
  std::vector<mps::Real> explicit_times;
  std::vector<mps::Real> implicit_times;
  std::vector<mps::Real> solve_times;

  void explicit_rhs(const mps::Real time_s, const std::span<const mps::Real> state,
                    const std::span<mps::Real> tendency) {
    explicit_times.push_back(time_s);
    tendency[0] = explicit_rate * state[0];
  }

  void implicit_rhs(const mps::Real time_s, const std::span<const mps::Real> state,
                    const std::span<mps::Real> tendency) {
    implicit_times.push_back(time_s);
    tendency[0] = implicit_rate * state[0];
  }

  void solve(const mps::Real time_s, const mps::Real implicit_time_s,
             const std::span<const mps::Real> right_hand_side,
             const std::span<mps::Real> result) {
    solve_times.push_back(time_s);
    result[0] = right_hand_side[0] / (1.0 - implicit_time_s * implicit_rate);
  }
};

[[nodiscard]] mps::Real integrate_split_scalar(const std::size_t steps) {
  mps::Ark2Imex stepper(1);
  SplitScalarProblem problem;
  problem.explicit_rate = 0.35;
  problem.implicit_rate = -1.25;
  std::array<mps::Real, 1> state{1.0};
  const mps::Real time_step_s = 1.0 / static_cast<mps::Real>(steps);
  for (std::size_t step = 0; step < steps; ++step) {
    stepper.step(
        static_cast<mps::Real>(step) * time_step_s, time_step_s, state,
        [&](const auto... arguments) { problem.explicit_rhs(arguments...); },
        [&](const auto... arguments) { problem.implicit_rhs(arguments...); },
        [&](const auto... arguments) { problem.solve(arguments...); });
  }
  return state[0];
}

}  // namespace

MPS_TEST_CASE("ARK2 coefficients satisfy the additive consistency conditions") {
  const auto table = mps::ark2_imex_table();
  MPS_CHECK_EQ(table.order, 2U);
  for (std::size_t stage = 0; stage < 3; ++stage) {
    mps::Real explicit_sum = 0.0;
    mps::Real implicit_sum = 0.0;
    for (std::size_t column = 0; column < 3; ++column) {
      explicit_sum += table.explicit_matrix[stage][column];
      implicit_sum += table.implicit_matrix[stage][column];
    }
    MPS_CHECK_NEAR(explicit_sum, table.nodes[stage], 2.0e-15);
    MPS_CHECK_NEAR(implicit_sum, table.nodes[stage], 2.0e-15);
  }
  mps::Real explicit_weight_sum = 0.0;
  mps::Real implicit_weight_sum = 0.0;
  for (std::size_t stage = 0; stage < 3; ++stage) {
    explicit_weight_sum += table.explicit_weights[stage];
    implicit_weight_sum += table.implicit_weights[stage];
  }
  MPS_CHECK_NEAR(explicit_weight_sum, 1.0, 2.0e-15);
  MPS_CHECK_NEAR(implicit_weight_sum, 1.0, 2.0e-15);
}

MPS_TEST_CASE("ARK2 uses the registered stage times and call counts") {
  mps::Ark2Imex stepper(1);
  SplitScalarProblem problem;
  problem.explicit_rate = 0.2;
  problem.implicit_rate = -0.4;
  std::array<mps::Real, 1> state{2.0};
  stepper.step(
      3.0, 0.5, state,
      [&](const auto... arguments) { problem.explicit_rhs(arguments...); },
      [&](const auto... arguments) { problem.implicit_rhs(arguments...); },
      [&](const auto... arguments) { problem.solve(arguments...); });

  const auto table = mps::ark2_imex_table();
  MPS_CHECK_EQ(problem.explicit_times.size(), 3U);
  MPS_CHECK_EQ(problem.implicit_times.size(), 3U);
  MPS_CHECK_EQ(problem.solve_times.size(), 2U);
  for (std::size_t stage = 0; stage < 3; ++stage) {
    const auto expected = 3.0 + 0.5 * table.nodes[stage];
    MPS_CHECK_NEAR(problem.explicit_times[stage], expected, 1.0e-15);
    MPS_CHECK_NEAR(problem.implicit_times[stage], expected, 1.0e-15);
    if (stage > 0) MPS_CHECK_NEAR(problem.solve_times[stage - 1], expected, 1.0e-15);
  }
}

MPS_TEST_CASE("ARK2 is second order for a split linear equation") {
  const auto exact = std::exp(-0.9);
  const auto coarse_error = std::abs(integrate_split_scalar(10) - exact);
  const auto medium_error = std::abs(integrate_split_scalar(20) - exact);
  const auto fine_error = std::abs(integrate_split_scalar(40) - exact);
  MPS_CHECK(coarse_error > medium_error);
  MPS_CHECK(medium_error > fine_error);
  MPS_CHECK(std::log2(coarse_error / medium_error) > 1.9);
  MPS_CHECK(std::log2(medium_error / fine_error) > 1.9);
}

MPS_TEST_CASE("ARK2 validates inputs and preserves state when a stage fails") {
  MPS_CHECK_THROWS_AS(mps::Ark2Imex(0), std::invalid_argument);
  mps::Ark2Imex stepper(1);
  std::array<mps::Real, 1> state{4.0};
  const auto zero = [](const mps::Real, const std::span<const mps::Real>,
                       const std::span<mps::Real> tendency) { tendency[0] = 0.0; };
  const auto failing_solve = [](const mps::Real, const mps::Real,
                                const std::span<const mps::Real>,
                                const std::span<mps::Real>) {
    throw std::runtime_error("comparison solve failed");
  };
  MPS_CHECK_THROWS_AS(stepper.step(0.0, 1.0, state, zero, zero, failing_solve),
                      std::runtime_error);
  MPS_CHECK_EQ(state[0], 4.0);
}

int main() { return mps::test::run_all(); }
