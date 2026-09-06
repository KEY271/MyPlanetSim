#include <array>
#include <atomic>
#include <cstdlib>
#include <limits>
#include <new>
#include <stdexcept>

#include "myplanetsim/numerics/gmres.hpp"
#include "support/test.hpp"

namespace {

std::atomic<bool> count_allocations{false};
std::atomic<std::size_t> allocation_count{0};

const mps::GmresLinearOperator matrix = [](const std::span<const mps::Real> input,
                                           const std::span<mps::Real> output) {
  output[0] = 4.0 * input[0] + input[1];
  output[1] = -2.0 * input[0] + 5.0 * input[1] + input[2];
  output[2] = input[1] + 3.0 * input[2];
};

}  // namespace

void* operator new(const std::size_t size) {
  if (count_allocations.load(std::memory_order_relaxed))
    allocation_count.fetch_add(1, std::memory_order_relaxed);
  if (void* pointer = std::malloc(size)) return pointer;
  throw std::bad_alloc();
}

void* operator new[](const std::size_t size) { return ::operator new(size); }
void operator delete(void* pointer) noexcept { std::free(pointer); }
void operator delete[](void* pointer) noexcept { std::free(pointer); }
void operator delete(void* pointer, std::size_t) noexcept { std::free(pointer); }
void operator delete[](void* pointer, std::size_t) noexcept { std::free(pointer); }

MPS_TEST_CASE("GMRES solves a nonsymmetric system") {
  const std::array<mps::Real, 3> expected{1.0, -2.0, 3.0};
  std::array<mps::Real, 3> right_hand_side{};
  matrix(expected, right_hand_side);
  std::array<mps::Real, 3> solution{};
  mps::GmresWorkspace workspace;
  const auto result = mps::restarted_gmres(matrix, right_hand_side, solution,
                                           {.restart = 3,
                                            .maximum_iterations = 6,
                                            .relative_tolerance = 1.0e-13,
                                            .absolute_tolerance = 1.0e-14},
                                           workspace);
  MPS_CHECK(result.converged());
  MPS_CHECK(result.iterations <= 3);
  MPS_CHECK(result.final_residual_norm <= 1.0e-12);
  for (std::size_t index = 0; index < solution.size(); ++index)
    MPS_CHECK_NEAR(solution[index], expected[index], 1.0e-12);
}

MPS_TEST_CASE("right preconditioning is applied to the solution basis") {
  const mps::GmresLinearOperator diagonal = [](const std::span<const mps::Real> input,
                                               const std::span<mps::Real> output) {
    output[0] = 1.0e-6 * input[0];
    output[1] = input[1];
    output[2] = 1.0e6 * input[2];
  };
  const mps::GmresPreconditioner inverse_diagonal =
      [](const std::span<const mps::Real> input, const std::span<mps::Real> output) {
        output[0] = 1.0e6 * input[0];
        output[1] = input[1];
        output[2] = 1.0e-6 * input[2];
      };
  const std::array<mps::Real, 3> right_hand_side{1.0, 2.0, 3.0};
  std::array<mps::Real, 3> solution{};
  mps::GmresWorkspace workspace;
  const auto result = mps::restarted_gmres(
      diagonal, right_hand_side, solution,
      {.restart = 1, .maximum_iterations = 1, .relative_tolerance = 1.0e-13}, workspace,
      inverse_diagonal);
  MPS_CHECK(result.converged());
  MPS_CHECK_EQ(result.iterations, 1U);
  MPS_CHECK_NEAR(solution[0], 1.0e6, 1.0e-7);
  MPS_CHECK_NEAR(solution[1], 2.0, 1.0e-12);
  MPS_CHECK_NEAR(solution[2], 3.0e-6, 1.0e-18);
}

MPS_TEST_CASE("zero residual converges without an iteration") {
  const mps::GmresLinearOperator identity = [](const std::span<const mps::Real> input,
                                               const std::span<mps::Real> output) {
    std::copy(input.begin(), input.end(), output.begin());
  };
  const std::array<mps::Real, 2> right_hand_side{2.0, -3.0};
  std::array<mps::Real, 2> solution{2.0, -3.0};
  mps::GmresWorkspace workspace;
  const auto result =
      mps::restarted_gmres(identity, right_hand_side, solution, {}, workspace);
  MPS_CHECK(result.converged());
  MPS_CHECK_EQ(result.iterations, 0U);
  MPS_CHECK_EQ(result.final_residual_norm, 0.0);
}

MPS_TEST_CASE("GMRES reports breakdown and iteration exhaustion") {
  const mps::GmresLinearOperator zero = [](const std::span<const mps::Real>,
                                           const std::span<mps::Real> output) {
    std::fill(output.begin(), output.end(), 0.0);
  };
  const std::array<mps::Real, 2> right_hand_side{1.0, 0.0};
  std::array<mps::Real, 2> solution{};
  mps::GmresWorkspace workspace;
  const auto breakdown =
      mps::restarted_gmres(zero, right_hand_side, solution, {}, workspace);
  MPS_CHECK(breakdown.status == mps::GmresStatus::kBreakdown);

  std::array<mps::Real, 3> exhausted_solution{};
  const auto exhausted = mps::restarted_gmres(
      matrix, std::array<mps::Real, 3>{6.0, -9.0, 7.0}, exhausted_solution,
      {.restart = 1,
       .maximum_iterations = 1,
       .relative_tolerance = 1.0e-15,
       .absolute_tolerance = 1.0e-15},
      workspace);
  MPS_CHECK(exhausted.status == mps::GmresStatus::kMaximumIterations);
}

MPS_TEST_CASE("GMRES rejects invalid and non-finite inputs") {
  std::array<mps::Real, 2> right_hand_side{1.0, 2.0};
  std::array<mps::Real, 2> solution{};
  mps::GmresWorkspace workspace;
  MPS_CHECK_THROWS_AS(
      mps::restarted_gmres({}, right_hand_side, solution, {}, workspace),
      std::invalid_argument);
  right_hand_side[1] = std::numeric_limits<mps::Real>::quiet_NaN();
  const mps::GmresLinearOperator identity = [](const std::span<const mps::Real> input,
                                               const std::span<mps::Real> output) {
    std::copy(input.begin(), input.end(), output.begin());
  };
  MPS_CHECK_THROWS_AS(
      mps::restarted_gmres(identity, right_hand_side, solution, {}, workspace),
      std::runtime_error);
}

MPS_TEST_CASE("warmed GMRES solve performs no heap allocations") {
  const std::array<mps::Real, 3> right_hand_side{6.0, -9.0, 7.0};
  std::array<mps::Real, 3> solution{};
  mps::GmresWorkspace workspace;
  const mps::GmresOptions options{
      .restart = 3, .maximum_iterations = 6, .relative_tolerance = 1.0e-12};
  static_cast<void>(
      mps::restarted_gmres(matrix, right_hand_side, solution, options, workspace));
  solution.fill(0.0);
  static_cast<void>(
      mps::restarted_gmres(matrix, right_hand_side, solution, options, workspace));
  solution.fill(0.0);

  allocation_count.store(0, std::memory_order_relaxed);
  count_allocations.store(true, std::memory_order_relaxed);
  const auto result =
      mps::restarted_gmres(matrix, right_hand_side, solution, options, workspace);
  count_allocations.store(false, std::memory_order_relaxed);

  MPS_CHECK(result.converged());
  MPS_CHECK_EQ(allocation_count.load(std::memory_order_relaxed), 0U);
}

int main() { return mps::test::run_all(); }
