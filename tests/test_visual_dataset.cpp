#include <chrono>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <string>

#include "myplanetsim/dynamics/dry_hydrostatic_driver.hpp"
#include "myplanetsim/io/run_metadata.hpp"
#include "myplanetsim/io/visual_dataset.hpp"
#include "support/test.hpp"

namespace {

mps::ExperimentConfig config() {
  return {.kind = mps::ExperimentKind::kDryHydrostatic,
          .planet = {2, 0, 10, 287, 1004, 100000},
          .run = {0, 1, .1, 0},
          .grid = {2},
          .vertical = {.levels = 2,
                       .a_half_pa = {1000, 500, 0},
                       .b_half = {0, .5, 1},
                       .surface_pressure_pa = 100000,
                       .minimum_surface_pressure_pa = 90000,
                       .maximum_surface_pressure_pa = 110000,
                       .minimum_pressure_thickness_pa = 100,
                       .initial_temperature_k = 280,
                       .initial_potential_temperature_k = 300,
                       .temperature_floor_k = 100,
                       .transport_scheme = mps::VerticalTransportScheme::kDonorCell,
                       .limiter = mps::VerticalLimiterKind::kNone,
                       .cfl = .5},
          .dry_hydrostatic = {},
          .diagnostics = {1},
          .statistics = {.enabled = true,
                         .start_time_s = 0,
                         .period_s = .5,
                         .fields = {"temperature", "pressure", "wind_speed"}},
          .output_directory = "unused"};
}

std::filesystem::path temporary_directory() {
  return std::filesystem::temp_directory_path() /
         ("myplanetsim-visual-" +
          std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
}

}  // namespace

MPS_TEST_CASE("visual dataset writes terrain, means, manifest, and restart sidecar") {
  const auto parameters = config();
  const mps::DryHydrostaticDriver driver(parameters);
  const auto state = driver.initial_state();
  const auto derived = driver.diagnose(state);
  const auto fields = mps::dry_visual_fields(parameters, driver.grid().cell_count(),
                                             derived.levels, false);
  const auto values =
      mps::dry_visual_values(parameters, driver.grid(), state, derived, fields);
  const auto directory = temporary_directory();
  std::filesystem::create_directories(directory);
  mps::VisualDatasetWriter writer(directory, parameters, driver.grid(),
                                  driver.orography(), {}, fields,
                                  mps::config_fingerprint(parameters));
  MPS_CHECK(std::filesystem::exists(directory / "terrain.bin"));
  MPS_CHECK(std::filesystem::exists(directory / "manifest.json"));
  writer.publish({.index = 0,
                  .scheduled_start_s = 0,
                  .scheduled_end_s = .5,
                  .actual_start_s = 0,
                  .actual_end_s = .25,
                  .weight_s = .25,
                  .complete = false,
                  .values = values});
  MPS_CHECK(std::filesystem::file_size(directory / "means/period_000000.bin") ==
            24 + 8 * values.size());
  std::ifstream manifest_input(directory / "manifest.json");
  const std::string manifest((std::istreambuf_iterator<char>(manifest_input)), {});
  MPS_CHECK(manifest.find("accepted_step_end_time_weighted_rectangle") !=
            std::string::npos);
  MPS_CHECK(manifest.find("\"complete\":false") != std::string::npos);

  const auto checkpoint = directory / "checkpoint.dat";
  {
    std::ofstream output(checkpoint);
    output << "checkpoint bytes";
  }
  mps::diagnostics::PeriodMeanAccumulator accumulator(values.size(), 0, .5);
  accumulator.observe(0, .25, values);
  const auto sidecar_path = mps::statistics_sidecar_path(checkpoint);
  mps::write_statistics_sidecar_file(sidecar_path,
                                     {.checkpoint_hash = mps::file_fnv1a64(checkpoint),
                                      .checkpoint_time_s = .25,
                                      .checkpoint_step = 1,
                                      .accumulator = accumulator.state(),
                                      .periods = writer.periods()});
  const auto restored = mps::read_statistics_sidecar_file(
      sidecar_path, mps::file_fnv1a64(checkpoint), .25, 1, values.size(), 0, .5);
  mps::diagnostics::PeriodMeanAccumulator resumed(restored.accumulator);
  MPS_CHECK_NEAR(resumed.partial()->values[0], values[0], 1e-14);
  MPS_CHECK_EQ(restored.periods.size(), 1U);
  std::filesystem::remove_all(directory);
}

int main() { return mps::test::run_all(); }
