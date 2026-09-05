#include <algorithm>
#include <filesystem>
#include <string>
#include <vector>

#include "myplanetsim/config/experiment_config.hpp"
#include "myplanetsim/dynamics/dry_hydrostatic_driver.hpp"
#include "support/test.hpp"

// Every shipped preset must at least start. `configs/phase6_dcmip_2_0_0.cfg` did not:
// its `a_half = [20544.8, 0, ..., 0]` with a uniform `b_half` put `p_half[1] = ps/15`
// below `p_half[0]`, so the coordinate failed validation and the Phase 6
// pressure-gradient measurement had to be taken on a substitute coordinate (ADR 0013).
// Nothing in CI noticed, because no gate ever loaded the preset.

namespace {

[[nodiscard]] std::vector<std::filesystem::path> preset_paths() {
  std::vector<std::filesystem::path> paths;
  for (const auto& entry : std::filesystem::directory_iterator(MPS_CONFIG_DIRECTORY)) {
    if (entry.path().extension() == ".cfg") paths.push_back(entry.path());
  }
  std::ranges::sort(paths);
  return paths;
}

}  // namespace

MPS_TEST_CASE("every shipped preset parses, validates, and starts") {
  const auto paths = preset_paths();
  MPS_CHECK(paths.size() >= 30);
  for (const auto& path : paths) {
    try {
      const auto config = mps::load_experiment_config(path.string());
      if (config.kind != mps::ExperimentKind::kDryHydrostatic) continue;
      const mps::DryHydrostaticDriver driver(config);
      const auto state = driver.initial_state();
      const auto rhs = driver.rhs(state);
      MPS_CHECK(rhs.horizontal_stable_time_step_s > 0.0);
    } catch (const std::exception& error) {
      ::mps::test::record_failure(__FILE__, __LINE__, path.filename().string(),
                                  error.what());
    }
  }
}

MPS_TEST_CASE("the DCMIP 2-0-0 coordinate is monotone") {
  const auto config = mps::load_experiment_config(
      (std::filesystem::path(MPS_CONFIG_DIRECTORY) / "phase6_dcmip_2_0_0.cfg")
          .string());
  const auto& a = config.vertical.a_half_pa;
  const auto& b = config.vertical.b_half;
  MPS_CHECK_EQ(a.size(), b.size());
  // Interfaces must increase downward at both the minimum and the maximum surface
  // pressure the preset admits, not only at its initial value.
  for (const mps::Real surface : {config.vertical.minimum_surface_pressure_pa,
                                  config.vertical.maximum_surface_pressure_pa}) {
    for (std::size_t k = 1; k < a.size(); ++k) {
      MPS_CHECK(a[k] + b[k] * surface > a[k - 1] + b[k - 1] * surface);
    }
  }
}

int main() { return mps::test::run_all(); }
