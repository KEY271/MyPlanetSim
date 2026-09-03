#include <sstream>
#include <stdexcept>
#include <string>

#include "myplanetsim/control/control_request.hpp"
#include "support/test.hpp"

namespace {

constexpr std::string_view kRequest = R"(
control.format_version = 1
control.run_id = run_01
control.end_time_s = 120
control.maximum_time_step_s = 2.5
control.frame_interval_steps = 4
control.frame_directory = runs/run_01/frames
initial_edits.count = 0
)";

constexpr std::string_view kEditedRequest = R"(
control.format_version = 1
control.run_id = run_02
control.end_time_s = 120
control.maximum_time_step_s = 2.5
control.frame_interval_steps = 4
control.frame_directory = runs/run_02/frames
initial_edits.count = 1
initial_edits.0.kind = gaussian_depth
initial_edits.0.center_x = 0
initial_edits.0.center_y = 0
initial_edits.0.center_z = 1
initial_edits.0.amplitude_m = 2
initial_edits.0.sigma_rad = 0.25
initial_edits.0.mass_policy = preserve_global
)";

[[nodiscard]] mps::ControlRequestV1 parse(const std::string_view text) {
  std::istringstream input{std::string(text)};
  return mps::parse_control_request(input);
}

}  // namespace

MPS_TEST_CASE("control request has a canonical round trip") {
  const auto first = parse(kRequest);
  std::ostringstream output;
  mps::write_control_request(output, first);
  const auto second = parse(output.str());
  MPS_CHECK_EQ(second.format_version, first.format_version);
  MPS_CHECK_EQ(second.run_id, first.run_id);
  MPS_CHECK_EQ(second.end_time_s, first.end_time_s);
  MPS_CHECK_EQ(second.maximum_time_step_s, first.maximum_time_step_s);
  MPS_CHECK_EQ(second.frame_interval_steps, first.frame_interval_steps);
  MPS_CHECK_EQ(second.frame_directory, first.frame_directory);
  MPS_CHECK_EQ(second.initial_edits.size(), first.initial_edits.size());
}

MPS_TEST_CASE("control request round trips a Gaussian edit") {
  const auto first = parse(kEditedRequest);
  MPS_CHECK_EQ(first.initial_edits.size(), std::size_t{1});
  MPS_CHECK_EQ(first.initial_edits[0].center_unit.z, 1.0);
  std::ostringstream output;
  mps::write_control_request(output, first);
  const auto second = parse(output.str());
  MPS_CHECK_EQ(second.initial_edits[0].amplitude_m, 2.0);
  MPS_CHECK(second.initial_edits[0].mass_policy == mps::MassPolicy::kPreserveGlobal);
}

MPS_TEST_CASE("control request rejects malformed structure") {
  MPS_CHECK_THROWS_AS(parse(std::string(kRequest) + "unknown.key = 1\n"),
                      std::runtime_error);
  MPS_CHECK_THROWS_AS(parse(std::string(kRequest) + "control.run_id = other\n"),
                      std::runtime_error);
  MPS_CHECK_THROWS_AS(parse("control.format_version = 1\n"), std::runtime_error);
}

MPS_TEST_CASE("control request rejects unsafe or out of range values") {
  const auto replace = [](std::string text, const std::string_view old_value,
                          const std::string_view new_value) {
    const auto position = text.find(old_value);
    text.replace(position, old_value.size(), new_value);
    return text;
  };
  auto invalid = std::string(kRequest);
  invalid = replace(invalid, "control.format_version = 1",
                    "control.format_version = 2");
  MPS_CHECK_THROWS_AS(parse(invalid), std::invalid_argument);

  invalid = replace(std::string(kRequest),
                    "control.frame_directory = runs/run_01/frames",
                    "control.frame_directory = ../escape");
  MPS_CHECK_THROWS_AS(parse(invalid), std::invalid_argument);

  invalid = replace(std::string(kRequest), "control.maximum_time_step_s = 2.5",
                    "control.maximum_time_step_s = nan");
  MPS_CHECK_THROWS_AS(parse(invalid), std::invalid_argument);
}

int main() { return mps::test::run_all(); }
