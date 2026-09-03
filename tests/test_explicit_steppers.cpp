#include <limits>
#include <span>
#include <stdexcept>
#include <vector>

#include "myplanetsim/core/types.hpp"
#include "myplanetsim/numerics/explicit_steppers.hpp"
#include "support/test.hpp"

namespace {

void zero_rhs(const mps::Real, const std::span<const mps::Real>,
              const std::span<mps::Real> tendency) {
  for (auto& value : tendency) {
    value = 0.0;
  }
}

void constant_rhs(const mps::Real, const std::span<const mps::Real>,
                  const std::span<mps::Real> tendency) {
  for (std::size_t index = 0; index < tendency.size(); ++index) {
    tendency[index] = static_cast<mps::Real>(index + 1);
  }
}

template <typename Stepper>
void check_common_stepper_behavior() {
  Stepper stepper(2);
  std::vector<mps::Real> state{1.0, -2.0};
  stepper.step(0.0, 0.25, state, zero_rhs);
  MPS_CHECK_EQ(state[0], 1.0);
  MPS_CHECK_EQ(state[1], -2.0);

  stepper.step(0.25, 0.5, state, constant_rhs);
  MPS_CHECK_NEAR(state[0], 1.5, 1.0e-15);
  MPS_CHECK_NEAR(state[1], -1.0, 1.0e-15);
}

}  // namespace

MPS_TEST_CASE("Forward Euler advances multiple components") {
  check_common_stepper_behavior<mps::ForwardEuler>();
}

MPS_TEST_CASE("SSP-RK3 advances multiple components") {
  check_common_stepper_behavior<mps::SspRk3>();
}

MPS_TEST_CASE("steppers reject invalid construction and arguments") {
  MPS_CHECK_THROWS_AS(mps::ForwardEuler(0), std::invalid_argument);
  MPS_CHECK_THROWS_AS(mps::SspRk3(0), std::invalid_argument);

  mps::ForwardEuler stepper(2);
  std::vector<mps::Real> state{1.0, 2.0};
  std::vector<mps::Real> wrong_size{1.0};
  MPS_CHECK_THROWS_AS(stepper.step(0.0, 0.0, state, zero_rhs), std::invalid_argument);
  MPS_CHECK_THROWS_AS(stepper.step(0.0, 1.0, wrong_size, zero_rhs),
                      std::invalid_argument);

  state[0] = std::numeric_limits<mps::Real>::infinity();
  MPS_CHECK_THROWS_AS(stepper.step(0.0, 1.0, state, zero_rhs), std::invalid_argument);
}

MPS_TEST_CASE("steppers reject non-finite tendencies before updating state") {
  mps::SspRk3 stepper(1);
  std::vector<mps::Real> state{3.0};
  const auto invalid_rhs = [](const mps::Real, const std::span<const mps::Real>,
                              const std::span<mps::Real> tendency) {
    tendency[0] = std::numeric_limits<mps::Real>::quiet_NaN();
  };

  MPS_CHECK_THROWS_AS(stepper.step(0.0, 0.1, state, invalid_rhs), std::runtime_error);
  MPS_CHECK_EQ(state[0], 3.0);
}

int main() { return mps::test::run_all(); }
