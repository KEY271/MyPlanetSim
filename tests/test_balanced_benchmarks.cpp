#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <string_view>

#include "myplanetsim/config/experiment_config.hpp"
#include "myplanetsim/dynamics/dry_hydrostatic_driver.hpp"
#include "support/test.hpp"

namespace {

struct ResidualNorms {
  mps::Real l2_m_s2;
  mps::Real linf_m_s2;
};

[[nodiscard]] ResidualNorms initial_momentum_residual(
    const std::string_view config_name, const mps::Index cells_per_panel) {
  auto config = mps::load_experiment_config(std::string(MPS_CONFIG_DIRECTORY) + "/" +
                                            std::string(config_name));
  config.grid.cells_per_panel = cells_per_panel;
  const mps::DryHydrostaticDriver driver(config);
  const auto state = driver.initial_state();
  const auto derived = driver.diagnose(state);
  const auto rhs = driver.rhs(state);

  mps::Real weighted_square = 0.0;
  mps::Real normalization = 0.0;
  mps::Real linf = 0.0;
  for (std::size_t cell = 0; cell < derived.cells; ++cell) {
    const mps::Real area = driver.grid().cells()[cell].area_m2;
    for (std::size_t level = 0; level < derived.levels; ++level) {
      const auto offset = mps::dry_hydrostatic_offset(cell, level, derived.levels);
      // M u is prognostic, so remove u dM/dt before dividing by M.
      const mps::Vec3 acceleration =
          (rhs.tendency.momentum[offset] -
           derived.velocity_m_s[offset] * rhs.tendency.air_mass[offset]) /
          derived.air_mass_kg_m2[offset];
      const mps::Real magnitude = mps::norm(acceleration);
      weighted_square += area * magnitude * magnitude;
      normalization += area;
      linf = std::max(linf, magnitude);
    }
  }
  return {.l2_m_s2 = std::sqrt(weighted_square / normalization), .linf_m_s2 = linf};
}

[[nodiscard]] std::array<ResidualNorms, 3> resolution_triplet(
    const std::string_view config_name) {
  std::array<ResidualNorms, 3> result{};
  for (std::size_t index = 0; index < result.size(); ++index) {
    result[index] = initial_momentum_residual(
        config_name, std::array<mps::Index, 3>{8, 16, 32}[index]);
  }
  return result;
}

void print_triplet(const std::string_view name,
                   const std::array<ResidualNorms, 3>& values) {
  std::cout << name;
  for (const auto& value : values) {
    std::cout << " l2=" << value.l2_m_s2 << " linf=" << value.linf_m_s2;
  }
  std::cout << '\n';
}

void require_decreasing(const std::array<ResidualNorms, 3>& values) {
  for (std::size_t index = 1; index < values.size(); ++index) {
    MPS_CHECK(values[index].l2_m_s2 < values[index - 1].l2_m_s2);
    MPS_CHECK(values[index].linf_m_s2 < values[index - 1].linf_m_s2);
  }
}

}  // namespace

MPS_TEST_CASE("UMJS14 steady initial residual decreases with resolution") {
  const auto values = resolution_triplet("phase5_umjs14_steady.cfg");
  print_triplet("umjs14", values);
  require_decreasing(values);
  // Registered after the Phase 9 analytic-state implementation. These are ceilings,
  // not claims of a pointwise exact discrete balance.
  MPS_CHECK(values[0].linf_m_s2 < 6.0e-4);
  MPS_CHECK(values[1].linf_m_s2 < 3.5e-4);
  MPS_CHECK(values[2].linf_m_s2 < 3.0e-4);
}

MPS_TEST_CASE("JW06 steady initial residual decreases with resolution") {
  const auto values = resolution_triplet("phase6_jw06_steady.cfg");
  print_triplet("jw06", values);
  require_decreasing(values);
  MPS_CHECK(values[0].linf_m_s2 < 8.0e-4);
  MPS_CHECK(values[1].linf_m_s2 < 4.5e-4);
  MPS_CHECK(values[2].linf_m_s2 < 3.3e-4);
}

MPS_TEST_CASE("DCMIP 2-0-0 pressure-gradient error is registered") {
  const auto values = resolution_triplet("phase6_dcmip_2_0_0.cfg");
  print_triplet("dcmip_2_0_0", values);
  // The repaired coordinate starts, but the current terrain-following pressure-gradient
  // discretization does not converge. ADR 0013 deliberately records this result for a
  // later discretization decision instead of weakening a convergence gate.
  MPS_CHECK(values[0].linf_m_s2 < 8.0e-4);
  MPS_CHECK(values[1].linf_m_s2 > 1.5e-3);
  MPS_CHECK(values[2].linf_m_s2 > 1.5e-3);
  MPS_CHECK(values[2].l2_m_s2 > values[1].l2_m_s2);
}

int main() { return mps::test::run_all(); }
