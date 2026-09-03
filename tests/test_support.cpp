#include "support/test.hpp"

#include <stdexcept>

MPS_TEST_CASE("boolean and equality checks") {
  MPS_CHECK(true);
  MPS_CHECK_EQ(2 + 2, 4);
  MPS_CHECK_NEAR(0.1 + 0.2, 0.3, 1.0e-12);
}

MPS_TEST_CASE("exception check") {
  MPS_CHECK_THROWS_AS(throw std::runtime_error("expected"), std::runtime_error);
}

int main() {
  return mps::test::run_all();
}
