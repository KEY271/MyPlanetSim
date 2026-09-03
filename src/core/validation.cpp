#include "myplanetsim/core/validation.hpp"

#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>

namespace mps {
namespace {

[[nodiscard]] std::invalid_argument invalid_value(const Real value,
                                                  const std::string_view parameter_name,
                                                  const std::string_view requirement) {
  std::ostringstream message;
  message << std::setprecision(std::numeric_limits<Real>::max_digits10)
          << parameter_name << "=" << value << " must be " << requirement;
  return std::invalid_argument(message.str());
}

}  // namespace

void require_finite(const Real value, const std::string_view parameter_name) {
  if (!std::isfinite(value)) {
    throw invalid_value(value, parameter_name, "finite");
  }
}

void require_positive(const Real value, const std::string_view parameter_name) {
  require_finite(value, parameter_name);
  if (value <= 0.0) {
    throw invalid_value(value, parameter_name, "greater than zero");
  }
}

void require_non_negative(const Real value, const std::string_view parameter_name) {
  require_finite(value, parameter_name);
  if (value < 0.0) {
    throw invalid_value(value, parameter_name, "non-negative");
  }
}

void require_in_closed_range(const Real value, const Real lower, const Real upper,
                             const std::string_view parameter_name) {
  require_finite(lower, "range lower bound");
  require_finite(upper, "range upper bound");
  if (lower > upper) {
    throw std::invalid_argument("range lower bound must not exceed upper bound");
  }
  require_finite(value, parameter_name);
  if (value < lower || value > upper) {
    std::ostringstream requirement;
    requirement << "in the closed range [" << lower << ", " << upper << ']';
    throw invalid_value(value, parameter_name, requirement.str());
  }
}

}  // namespace mps
