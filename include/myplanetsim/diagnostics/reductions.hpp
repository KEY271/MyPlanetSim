#pragma once

#include <span>

#include "myplanetsim/core/types.hpp"

namespace mps::diagnostics {

struct MinMax {
  Real minimum;
  Real maximum;
};

struct ErrorNorms {
  Real l1;
  Real l2;
  Real linf;
};

class CompensatedAccumulator {
 public:
  void add(Real value);
  [[nodiscard]] Real value() const;

 private:
  Real sum_ = 0.0;
  Real correction_ = 0.0;
};

[[nodiscard]] bool all_finite(std::span<const Real> values) noexcept;
[[nodiscard]] Real compensated_sum(std::span<const Real> values);
[[nodiscard]] MinMax min_max(std::span<const Real> values);
[[nodiscard]] ErrorNorms weighted_error_norms(std::span<const Real> actual,
                                              std::span<const Real> expected,
                                              std::span<const Real> weights);
[[nodiscard]] Real relative_drift(Real current, Real initial, Real scale);

}  // namespace mps::diagnostics
