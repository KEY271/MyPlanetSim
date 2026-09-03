#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <sstream>

#include "myplanetsim/config/experiment_config.hpp"
#include "myplanetsim/io/checkpoint.hpp"
#include "myplanetsim/io/run_metadata.hpp"
#include "myplanetsim/phase0/ode_experiment.hpp"
#include "support/test.hpp"

namespace {

[[nodiscard]] mps::ExperimentConfig config_with_time_step(const mps::Real time_step_s) {
  return mps::ExperimentConfig{
      .planet = mps::PlanetParameters::earth_like(),
      .run = {.start_time_s = 0.0,
              .end_time_s = 1.0,
              .time_step_s = time_step_s,
              .random_seed = 5489},
      .ode = {.initial_value = 1.0, .decay_rate_s_1 = 1.0},
      .output_directory = "output",
  };
}

[[nodiscard]] mps::Real observed_order(const mps::Real coarse_error,
                                       const mps::Real fine_error) {
  return std::log(coarse_error / fine_error) / std::log(2.0);
}

template <std::size_t Size>
void check_convergence(const mps::IntegratorKind integrator,
                       const std::array<mps::Real, Size>& time_steps,
                       const mps::Real minimum_order) {
  std::array<mps::Real, Size> errors{};
  for (std::size_t index = 0; index < Size; ++index) {
    errors[index] =
        mps::run_ode_experiment(config_with_time_step(time_steps[index]), integrator)
            .absolute_error;
  }

  std::cout << mps::integrator_name(integrator) << " convergence:\n";
  for (std::size_t index = 0; index < Size; ++index) {
    std::cout << "  dt=" << time_steps[index] << " error=" << errors[index];
    if (index > 0) {
      std::cout << " order=" << observed_order(errors[index - 1], errors[index]);
    }
    std::cout << '\n';
  }

  for (std::size_t index = 1; index < Size; ++index) {
    MPS_CHECK(errors[index] < errors[index - 1]);
    MPS_CHECK(observed_order(errors[index - 1], errors[index]) >= minimum_order);
  }
}

}  // namespace

MPS_TEST_CASE("Forward Euler and SSP-RK3 meet temporal convergence gates") {
  constexpr std::array<mps::Real, 4> time_steps{0.2, 0.1, 0.05, 0.025};
  check_convergence(mps::IntegratorKind::kForwardEuler, time_steps, 0.7);
  check_convergence(mps::IntegratorKind::kSspRk3, time_steps, 2.7);
}

MPS_TEST_CASE("serialized restart matches continuous integration bitwise") {
  const auto config = config_with_time_step(0.1);
  const auto fingerprint = mps::config_fingerprint(config);
  const auto continuous = mps::run_ode_experiment(config, mps::IntegratorKind::kSspRk3);
  const auto partial =
      mps::run_ode_experiment(config, mps::IntegratorKind::kSspRk3, std::nullopt, 5);

  std::ostringstream serialized;
  mps::write_checkpoint(serialized, mps::Checkpoint{.time_s = partial.state.time_s,
                                                    .step = partial.state.step,
                                                    .state = {partial.state.value},
                                                    .config_fingerprint = fingerprint});
  std::istringstream input(serialized.str());
  const auto checkpoint = mps::read_checkpoint(input, fingerprint);
  const auto restarted =
      mps::run_ode_experiment(config, mps::IntegratorKind::kSspRk3,
                              mps::OdeState{.time_s = checkpoint.time_s,
                                            .step = checkpoint.step,
                                            .value = checkpoint.state.front()});

  MPS_CHECK_EQ(restarted.state.time_s, continuous.state.time_s);
  MPS_CHECK_EQ(restarted.state.step, continuous.state.step);
  MPS_CHECK_EQ(std::bit_cast<std::uint64_t>(restarted.state.value),
               std::bit_cast<std::uint64_t>(continuous.state.value));
}

MPS_TEST_CASE("configuration and seed retain deterministic identities") {
  const auto config = config_with_time_step(0.1);
  MPS_CHECK_EQ(mps::canonical_config_text(config), mps::canonical_config_text(config));
  MPS_CHECK_EQ(mps::config_fingerprint(config), mps::config_fingerprint(config));

  auto first = mps::make_random_engine(config.run.random_seed);
  auto second = mps::make_random_engine(config.run.random_seed);
  for (int draw = 0; draw < 32; ++draw) {
    MPS_CHECK_EQ(first(), second());
  }
}

int main() { return mps::test::run_all(); }
