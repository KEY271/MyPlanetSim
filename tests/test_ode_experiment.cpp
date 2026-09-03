#include "myplanetsim/config/experiment_config.hpp"
#include "myplanetsim/phase0/ode_experiment.hpp"
#include "support/test.hpp"

namespace {

[[nodiscard]] mps::ExperimentConfig sample_config() {
  return mps::ExperimentConfig{
      .planet = mps::PlanetParameters::earth_like(),
      .run = {.start_time_s = 0.0,
              .end_time_s = 1.0,
              .time_step_s = 0.3,
              .random_seed = 1},
      .ode = {.initial_value = 1.0, .decay_rate_s_1 = 1.0},
      .output_directory = "output",
  };
}

}  // namespace

MPS_TEST_CASE("ODE experiment lands exactly on its final time") {
  const auto config = sample_config();
  const auto result =
      mps::run_ode_experiment(config, mps::IntegratorKind::kForwardEuler);
  MPS_CHECK(result.reached_end_time);
  MPS_CHECK_EQ(result.state.time_s, config.run.end_time_s);
  MPS_CHECK_EQ(result.state.step, 4U);
  MPS_CHECK(result.absolute_error > 0.0);
}

MPS_TEST_CASE("ODE experiment does not add a roundoff-sized final step") {
  auto config = sample_config();
  config.run.time_step_s = 0.1;
  const auto result = mps::run_ode_experiment(config, mps::IntegratorKind::kSspRk3);
  MPS_CHECK_EQ(result.state.time_s, 1.0);
  MPS_CHECK_EQ(result.state.step, 10U);
}

MPS_TEST_CASE("integrator names have a strict round trip") {
  MPS_CHECK(mps::parse_integrator("euler") == mps::IntegratorKind::kForwardEuler);
  MPS_CHECK(mps::parse_integrator("ssprk3") == mps::IntegratorKind::kSspRk3);
  MPS_CHECK_EQ(mps::integrator_name(mps::IntegratorKind::kSspRk3), "ssprk3");
  MPS_CHECK_THROWS_AS(mps::parse_integrator("rk4"), std::invalid_argument);
}

int main() { return mps::test::run_all(); }
