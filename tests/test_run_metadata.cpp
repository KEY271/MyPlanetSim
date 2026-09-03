#include <sstream>
#include <string>

#include "myplanetsim/config/experiment_config.hpp"
#include "myplanetsim/io/run_metadata.hpp"
#include "support/test.hpp"

namespace {

[[nodiscard]] mps::ExperimentConfig sample_config() {
  return mps::ExperimentConfig{
      .planet = mps::PlanetParameters::earth_like(),
      .run =
          mps::RunParameters{
              .start_time_s = 0.0,
              .end_time_s = 1.0,
              .time_step_s = 0.1,
              .random_seed = 1234,
          },
      .ode =
          mps::OdeParameters{
              .initial_value = 1.0,
              .decay_rate_s_1 = 1.0,
          },
      .output_directory = "output",
  };
}

}  // namespace

MPS_TEST_CASE("configuration fingerprint is deterministic and sensitive") {
  auto config = sample_config();
  const auto first = mps::config_fingerprint(config);
  MPS_CHECK_EQ(first, mps::config_fingerprint(config));
  MPS_CHECK_EQ(first.size(), 16U);

  config.run.time_step_s = 0.05;
  MPS_CHECK(first != mps::config_fingerprint(config));
}

MPS_TEST_CASE("run metadata contains required canonical fields") {
  const auto config = sample_config();
  const auto metadata = mps::make_run_metadata(config);
  std::ostringstream output;
  mps::write_run_metadata(output, metadata, config);
  const auto text = output.str();

  MPS_CHECK(text.find("metadata.format_version = 1\n") != std::string::npos);
  MPS_CHECK(text.find("metadata.model_version = 0.1.0\n") != std::string::npos);
  MPS_CHECK(text.find("metadata.compiler_id = ") != std::string::npos);
  MPS_CHECK(text.find("metadata.random_seed = 1234\n") != std::string::npos);
  MPS_CHECK(text.find("metadata.config_fingerprint = " +
                      mps::config_fingerprint(config) + "\n") != std::string::npos);
  MPS_CHECK(text.find("configuration.begin\n") != std::string::npos);
  MPS_CHECK(text.ends_with("configuration.end\n"));
}

MPS_TEST_CASE("explicit seeds produce reproducible random sequences") {
  auto first = mps::make_random_engine(77);
  auto second = mps::make_random_engine(77);
  auto different = mps::make_random_engine(78);

  for (int draw = 0; draw < 8; ++draw) {
    MPS_CHECK_EQ(first(), second());
  }
  MPS_CHECK(first() != different());
}

int main() { return mps::test::run_all(); }
