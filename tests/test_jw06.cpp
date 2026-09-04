#include <algorithm>
#include <cmath>

#include "myplanetsim/dynamics/dry_hydrostatic_driver.hpp"
#include "support/test.hpp"

namespace {
mps::ExperimentConfig config(mps::DryHydrostaticTestCase test_case) {
  return {.kind = mps::ExperimentKind::kDryHydrostatic,
          .planet = {6.371229e6, 7.29212e-5, 9.80616, 287, 1004.5, 100000},
          .run = {0, 60, 30, 0},
          .grid = {8},
          .vertical = {.levels = 4,
                       .a_half_pa = {1000, 0, 0, 0, 0},
                       .b_half = {0, .25, .5, .75, 1},
                       .surface_pressure_pa = 100000,
                       .minimum_surface_pressure_pa = 99000,
                       .maximum_surface_pressure_pa = 101000,
                       .minimum_pressure_thickness_pa = 100,
                       .initial_temperature_k = 288,
                       .initial_potential_temperature_k = 300,
                       .temperature_floor_k = 100,
                       .transport_scheme = mps::VerticalTransportScheme::kLinear,
                       .limiter = mps::VerticalLimiterKind::kMinmod,
                       .cfl = .5},
          .dry_hydrostatic = {.test_case = test_case},
          .orography = {.kind = mps::OrographyKind::kJw06},
          .diagnostics = {1},
          .output_directory = "x"};
}
}  // namespace

MPS_TEST_CASE("JW06 steady and perturbed initial states form a paired regression") {
  const mps::DryHydrostaticDriver steady_driver(
      config(mps::DryHydrostaticTestCase::kJw06Steady));
  const mps::DryHydrostaticDriver perturbed_driver(
      config(mps::DryHydrostaticTestCase::kJw06Baroclinic));
  const auto steady = steady_driver.initial_state();
  const auto perturbed = perturbed_driver.initial_state();
  MPS_CHECK(steady.surface_pressure_pa == perturbed.surface_pressure_pa);
  double maximum_difference = 0;
  for (std::size_t offset = 0; offset < steady.horizontal_momentum_mass_kg_m_s.size();
       ++offset) {
    maximum_difference =
        std::max(maximum_difference,
                 mps::norm(perturbed.horizontal_momentum_mass_kg_m_s[offset] -
                           steady.horizontal_momentum_mass_kg_m_s[offset]));
  }
  MPS_CHECK(maximum_difference > 0);
  const auto derived = steady_driver.diagnose(steady);
  for (const auto temperature : derived.temperature_k) {
    MPS_CHECK(std::isfinite(temperature));
    MPS_CHECK(temperature > 180);
    MPS_CHECK(temperature < 330);
  }
}

int main() { return mps::test::run_all(); }
