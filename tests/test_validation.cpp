#include <limits>
#include <stdexcept>
#include <string>

#include "myplanetsim/core/types.hpp"
#include "myplanetsim/core/validation.hpp"
#include "support/test.hpp"

MPS_TEST_CASE("valid values pass validation") {
  mps::require_finite(0.0, "zero");
  mps::require_positive(1.0, "positive");
  mps::require_non_negative(0.0, "non_negative");
  mps::require_in_closed_range(-1.0, -1.0, 1.0, "lower_boundary");
  mps::require_in_closed_range(1.0, -1.0, 1.0, "upper_boundary");
}

MPS_TEST_CASE("invalid finite and sign values are rejected") {
  const auto infinity = std::numeric_limits<mps::Real>::infinity();
  const auto nan = std::numeric_limits<mps::Real>::quiet_NaN();

  MPS_CHECK_THROWS_AS(mps::require_finite(infinity, "infinity"), std::invalid_argument);
  MPS_CHECK_THROWS_AS(mps::require_finite(nan, "nan"), std::invalid_argument);
  MPS_CHECK_THROWS_AS(mps::require_positive(0.0, "zero"), std::invalid_argument);
  MPS_CHECK_THROWS_AS(mps::require_non_negative(-1.0, "negative"),
                      std::invalid_argument);
  MPS_CHECK_THROWS_AS(mps::require_in_closed_range(2.0, -1.0, 1.0, "outside"),
                      std::invalid_argument);
  MPS_CHECK_THROWS_AS(mps::require_in_closed_range(0.0, 1.0, -1.0, "bad_range"),
                      std::invalid_argument);
}

MPS_TEST_CASE("validation errors identify the parameter") {
  try {
    mps::require_positive(-2.0, "planet.radius_m");
    MPS_CHECK(false);
  } catch (const std::invalid_argument& error) {
    MPS_CHECK(std::string(error.what()).find("planet.radius_m") != std::string::npos);
    MPS_CHECK(std::string(error.what()).find("-2") != std::string::npos);
  }
}

int main() { return mps::test::run_all(); }
