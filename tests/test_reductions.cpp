#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "myplanetsim/diagnostics/reductions.hpp"
#include "support/test.hpp"

MPS_TEST_CASE("compensated sum recovers a small residual") {
  const std::array<mps::Real, 3> values{1.0e16, 1.0, -1.0e16};
  MPS_CHECK_EQ(mps::diagnostics::compensated_sum(values), 1.0);
  MPS_CHECK_EQ(mps::diagnostics::compensated_sum(values),
               mps::diagnostics::compensated_sum(values));
}

MPS_TEST_CASE("min max and finite checks report field bounds") {
  const std::array<mps::Real, 4> values{4.0, -2.0, 7.0, 1.0};
  const auto bounds = mps::diagnostics::min_max(values);
  MPS_CHECK_EQ(bounds.minimum, -2.0);
  MPS_CHECK_EQ(bounds.maximum, 7.0);
  MPS_CHECK(mps::diagnostics::all_finite(values));

  const std::array<mps::Real, 1> invalid{std::numeric_limits<mps::Real>::infinity()};
  MPS_CHECK(!mps::diagnostics::all_finite(invalid));
  MPS_CHECK_THROWS_AS(mps::diagnostics::min_max(invalid), std::invalid_argument);
}

MPS_TEST_CASE("weighted error norms match hand calculations") {
  const std::array<mps::Real, 2> actual{2.0, 4.0};
  const std::array<mps::Real, 2> expected{1.0, 2.0};
  const std::array<mps::Real, 2> weights{1.0, 3.0};
  const auto norms = mps::diagnostics::weighted_error_norms(actual, expected, weights);

  MPS_CHECK_NEAR(norms.l1, 1.75, 1.0e-15);
  MPS_CHECK_NEAR(norms.l2, std::sqrt(13.0 / 4.0), 1.0e-15);
  MPS_CHECK_EQ(norms.linf, 2.0);
}

MPS_TEST_CASE("diagnostics reject invalid inputs") {
  const std::array<mps::Real, 1> one{1.0};
  const std::array<mps::Real, 2> two{1.0, 2.0};
  const std::array<mps::Real, 1> negative_weight{-1.0};
  const std::array<mps::Real, 0> empty{};

  MPS_CHECK_THROWS_AS(mps::diagnostics::min_max(empty), std::invalid_argument);
  MPS_CHECK_THROWS_AS(mps::diagnostics::weighted_error_norms(one, two, one),
                      std::invalid_argument);
  MPS_CHECK_THROWS_AS(mps::diagnostics::weighted_error_norms(one, one, negative_weight),
                      std::invalid_argument);
  MPS_CHECK_THROWS_AS(mps::diagnostics::relative_drift(1.0, 1.0, 0.0),
                      std::invalid_argument);
}

MPS_TEST_CASE("relative drift handles a zero initial value") {
  MPS_CHECK_NEAR(mps::diagnostics::relative_drift(0.25, 0.0, 0.5), 0.5, 1.0e-15);
  MPS_CHECK_NEAR(mps::diagnostics::relative_drift(11.0, 10.0, 1.0), 0.1, 1.0e-15);
}

int main() { return mps::test::run_all(); }
