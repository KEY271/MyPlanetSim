#include <bit>
#include <cstdint>
#include <filesystem>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>

#include "myplanetsim/io/checkpoint.hpp"
#include "support/test.hpp"

namespace {

[[nodiscard]] mps::Checkpoint sample_checkpoint() {
  return mps::Checkpoint{
      .time_s = 0.75,
      .step = 7,
      .state = {0.1, -0.0, std::numeric_limits<mps::Real>::min(),
                std::numeric_limits<mps::Real>::max()},
      .config_fingerprint = "0123456789abcdef",
  };
}

}  // namespace

MPS_TEST_CASE("checkpoint round trips every binary64 bit pattern used") {
  const auto original = sample_checkpoint();
  std::ostringstream output;
  mps::write_checkpoint(output, original);
  std::istringstream input(output.str());
  const auto restored = mps::read_checkpoint(input, original.config_fingerprint);

  MPS_CHECK_EQ(restored.time_s, original.time_s);
  MPS_CHECK_EQ(restored.step, original.step);
  MPS_CHECK_EQ(restored.state.size(), original.state.size());
  for (std::size_t index = 0; index < original.state.size(); ++index) {
    MPS_CHECK_EQ(std::bit_cast<std::uint64_t>(restored.state[index]),
                 std::bit_cast<std::uint64_t>(original.state[index]));
  }
}

MPS_TEST_CASE("checkpoint rejects mismatches corruption and truncation") {
  const auto checkpoint = sample_checkpoint();
  std::ostringstream output;
  mps::write_checkpoint(output, checkpoint);

  std::istringstream wrong_config(output.str());
  MPS_CHECK_THROWS_AS(mps::read_checkpoint(wrong_config, "different"),
                      std::runtime_error);

  auto wrong_version_text = output.str();
  wrong_version_text.replace(wrong_version_text.find("version = 1"), 11, "version = 2");
  std::istringstream wrong_version(wrong_version_text);
  MPS_CHECK_THROWS_AS(
      mps::read_checkpoint(wrong_version, checkpoint.config_fingerprint),
      std::runtime_error);

  auto truncated_text = output.str();
  truncated_text.resize(truncated_text.find("state.end"));
  std::istringstream truncated(truncated_text);
  MPS_CHECK_THROWS_AS(mps::read_checkpoint(truncated, checkpoint.config_fingerprint),
                      std::runtime_error);
}

MPS_TEST_CASE("checkpoint file is promoted from a temporary file") {
  const auto path =
      std::filesystem::temp_directory_path() / "myplanetsim_phase0_checkpoint_test.txt";
  const auto temporary_path = std::filesystem::path(path.string() + ".tmp");
  std::error_code ignored;
  std::filesystem::remove(path, ignored);
  std::filesystem::remove(temporary_path, ignored);

  const auto original = sample_checkpoint();
  mps::write_checkpoint_file(path, original);
  MPS_CHECK(std::filesystem::exists(path));
  MPS_CHECK(!std::filesystem::exists(temporary_path));
  const auto restored = mps::read_checkpoint_file(path, original.config_fingerprint);
  MPS_CHECK_EQ(restored.step, original.step);

  std::filesystem::remove(path, ignored);
}

int main() { return mps::test::run_all(); }
