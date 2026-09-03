#include "myplanetsim/diagnostics/reductions.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "myplanetsim/core/validation.hpp"

namespace mps::diagnostics {

bool all_finite(const std::span<const Real> values) noexcept {
  return std::ranges::all_of(values,
                             [](const Real value) { return std::isfinite(value); });
}

Real compensated_sum(const std::span<const Real> values) {
  Real sum = 0.0;
  Real correction = 0.0;
  for (const Real value : values) {
    require_finite(value, "reduction value");
    const Real next = sum + value;
    if (std::abs(sum) >= std::abs(value)) {
      correction += (sum - next) + value;
    } else {
      correction += (value - next) + sum;
    }
    sum = next;
  }
  const Real result = sum + correction;
  require_finite(result, "reduction result");
  return result;
}

MinMax min_max(const std::span<const Real> values) {
  if (values.empty()) {
    throw std::invalid_argument("min_max requires at least one value");
  }
  if (!all_finite(values)) {
    throw std::invalid_argument("min_max input contains a non-finite value");
  }
  const auto [minimum, maximum] = std::ranges::minmax_element(values);
  return MinMax{.minimum = *minimum, .maximum = *maximum};
}

ErrorNorms weighted_error_norms(const std::span<const Real> actual,
                                const std::span<const Real> expected,
                                const std::span<const Real> weights) {
  if (actual.empty()) {
    throw std::invalid_argument("error norms require at least one value");
  }
  if (actual.size() != expected.size() || actual.size() != weights.size()) {
    throw std::invalid_argument("error norm inputs must have equal sizes");
  }

  std::vector<Real> weighted_absolute_error(actual.size());
  std::vector<Real> weighted_squared_error(actual.size());
  Real maximum_error = 0.0;
  for (std::size_t index = 0; index < actual.size(); ++index) {
    require_finite(actual[index], "actual value");
    require_finite(expected[index], "expected value");
    require_positive(weights[index], "error norm weight");
    const Real error = actual[index] - expected[index];
    const Real absolute_error = std::abs(error);
    weighted_absolute_error[index] = weights[index] * absolute_error;
    weighted_squared_error[index] = weights[index] * error * error;
    maximum_error = std::max(maximum_error, absolute_error);
  }

  const Real weight_sum = compensated_sum(weights);
  const Real l1 = compensated_sum(weighted_absolute_error) / weight_sum;
  const Real l2 = std::sqrt(compensated_sum(weighted_squared_error) / weight_sum);
  require_finite(l1, "L1 error norm");
  require_finite(l2, "L2 error norm");
  return ErrorNorms{.l1 = l1, .l2 = l2, .linf = maximum_error};
}

Real relative_drift(const Real current, const Real initial, const Real scale) {
  require_finite(current, "current conserved quantity");
  require_finite(initial, "initial conserved quantity");
  require_positive(scale, "conservation scale");
  return (current - initial) / std::max(std::abs(initial), scale);
}

}  // namespace mps::diagnostics
