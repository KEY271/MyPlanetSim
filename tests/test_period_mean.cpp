#include <array>
#include <limits>
#include <stdexcept>

#include "myplanetsim/diagnostics/period_mean.hpp"
#include "support/test.hpp"

MPS_TEST_CASE("period means split variable steps across half-open windows") {
  mps::diagnostics::PeriodMeanAccumulator mean(2, 2.0, 4.0);
  mean.observe(0.0, 3.0, std::array{2.0, 20.0});
  mean.observe(3.0, 8.0, std::array{4.0, 40.0});
  const auto completed = mean.take_completed();
  MPS_CHECK(completed.size() == 1);
  MPS_CHECK(completed[0].index == 0);
  MPS_CHECK_NEAR(completed[0].weight_s, 4.0, 1e-14);
  MPS_CHECK_NEAR(completed[0].values[0], 3.5, 1e-14);
  MPS_CHECK_NEAR(completed[0].values[1], 35.0, 1e-14);
  MPS_CHECK(completed[0].complete);
  const auto partial = mean.partial();
  MPS_CHECK(partial.has_value());
  MPS_CHECK(partial->index == 1);
  MPS_CHECK_NEAR(partial->weight_s, 2.0, 1e-14);
  MPS_CHECK_NEAR(partial->values[0], 4.0, 1e-14);
}

MPS_TEST_CASE("period means handle one accepted step spanning several windows") {
  mps::diagnostics::PeriodMeanAccumulator mean(1, 1.0, 2.0);
  mean.observe(0.0, 8.0, std::array{7.0});
  const auto completed = mean.take_completed();
  MPS_CHECK(completed.size() == 3);
  for (const auto& window : completed) {
    MPS_CHECK_NEAR(window.weight_s, 2.0, 1e-14);
    MPS_CHECK_NEAR(window.values[0], 7.0, 1e-14);
  }
  const auto partial = mean.partial();
  MPS_CHECK(partial.has_value());
  MPS_CHECK(partial->index == 3);
  MPS_CHECK_NEAR(partial->weight_s, 1.0, 1e-14);
}

MPS_TEST_CASE("period means reject bad accepted-step notifications") {
  mps::diagnostics::PeriodMeanAccumulator mean(1, 0.0, 2.0);
  mean.observe(0.0, 1.0, std::array{1.0});
  MPS_CHECK_THROWS_AS(mean.observe(0.5, 1.5, std::array{1.0}), std::invalid_argument);
  MPS_CHECK_THROWS_AS(mean.observe(2.0, 1.0, std::array{1.0}), std::invalid_argument);
  MPS_CHECK_THROWS_AS(
      mean.observe(1.0, 2.0, std::array{std::numeric_limits<double>::quiet_NaN()}),
      std::invalid_argument);
}

MPS_TEST_CASE("zero-duration initial and duplicate notifications do not contribute") {
  mps::diagnostics::PeriodMeanAccumulator mean(1, 0.0, 2.0);
  mean.observe(0.0, 0.0, std::array{99.0});
  mean.observe(0.0, 1.0, std::array{3.0});
  mean.observe(1.0, 1.0, std::array{99.0});
  const auto partial = mean.partial();
  MPS_CHECK(partial.has_value());
  MPS_CHECK_NEAR(partial->weight_s, 1.0, 1e-14);
  MPS_CHECK_NEAR(partial->values[0], 3.0, 1e-14);
}

int main() { return mps::test::run_all(); }
