#include <algorithm>
#include <cmath>
#include <cstddef>
#include <random>
#include <stdexcept>
#include <vector>

#include "myplanetsim/physics/dry_convective_adjustment.hpp"
#include "myplanetsim/vertical/hydrostatic_column.hpp"
#include "support/test.hpp"

namespace {

[[nodiscard]] mps::DryConvectiveAdjustmentInput input(
    const std::vector<mps::Real>& theta, const std::vector<mps::Real>& mass,
    const std::vector<mps::Real>& exner_full,
    const std::vector<mps::Real>& exner_half) {
  return {.potential_temperature_k = theta,
          .air_mass_kg_m2 = mass,
          .exner_full = exner_full,
          .exner_half = exner_half,
          .heat_capacity_cp_j_kg_k = 1004.0,
          .heat_capacity_cv_j_kg_k = 717.0,
          .stability_tolerance_k = 1e-10};
}

[[nodiscard]] mps::Real enthalpy(const std::vector<mps::Real>& theta,
                                 const std::vector<mps::Real>& mass,
                                 const std::vector<mps::Real>& exner) {
  mps::Real result = 0.0;
  for (std::size_t level = 0; level < theta.size(); ++level)
    result += 1004.0 * mass[level] * exner[level] * theta[level];
  return result;
}

}  // namespace

MPS_TEST_CASE("two-layer adjustment has the enthalpy-weighted analytic solution") {
  const std::vector<mps::Real> theta{280.0, 320.0};
  const std::vector<mps::Real> mass{2.0, 1.0};
  const std::vector<mps::Real> exner{0.5, 1.0};
  const std::vector<mps::Real> exner_half{0.4, 0.7, 1.0};
  const auto result =
      mps::dry_convective_adjustment(input(theta, mass, exner, exner_half));
  MPS_CHECK_NEAR(result.adjusted_potential_temperature_k[0], 300.0, 1e-13);
  MPS_CHECK_NEAR(result.adjusted_potential_temperature_k[1], 300.0, 1e-13);
  MPS_CHECK_EQ(result.diagnostics.adjusted_block_count, 1U);
  MPS_CHECK_EQ(result.diagnostics.adjusted_layer_count, 2U);
  MPS_CHECK_NEAR(result.diagnostics.enthalpy_change_j_m2, 0.0, 1e-8);
  MPS_CHECK(result.diagnostics.unstable_interface_fraction_before == 1.0);
  MPS_CHECK(result.diagnostics.unstable_interface_fraction_after == 0.0);
}

MPS_TEST_CASE("PAVA revisits upper blocks and preserves disconnected regions") {
  const std::vector<mps::Real> theta{330.0, 300.0, 310.0, 280.0, 290.0};
  const std::vector<mps::Real> mass(theta.size(), 1.0);
  const std::vector<mps::Real> exner(theta.size(), 1.0);
  const std::vector<mps::Real> exner_half{0.2, 0.35, 0.5, 0.65, 0.8, 1.0};
  const auto result =
      mps::dry_convective_adjustment(input(theta, mass, exner, exner_half));
  const std::vector<mps::Real> expected{330.0, 305.0, 305.0, 285.0, 285.0};
  for (std::size_t level = 0; level < theta.size(); ++level)
    MPS_CHECK_NEAR(result.adjusted_potential_temperature_k[level], expected[level],
                   1e-13);
  MPS_CHECK_EQ(result.diagnostics.adjusted_block_count, 2U);

  const std::vector<mps::Real> cascading{300.0, 280.0, 340.0};
  const std::vector<mps::Real> cascading_mass(cascading.size(), 1.0);
  const std::vector<mps::Real> cascading_exner(cascading.size(), 1.0);
  const std::vector<mps::Real> cascading_half{0.3, 0.5, 0.7, 1.0};
  const auto merged = mps::dry_convective_adjustment(
      input(cascading, cascading_mass, cascading_exner, cascading_half));
  for (const auto value : merged.adjusted_potential_temperature_k)
    MPS_CHECK_NEAR(value, 920.0 / 3.0, 1e-13);
}

MPS_TEST_CASE("stable columns are identities and adjustment is idempotent") {
  const std::vector<mps::Real> stable{330.0, 320.0, 310.0, 300.0};
  const std::vector<mps::Real> mass{1.0, 2.0, 3.0, 4.0};
  const std::vector<mps::Real> exner{0.5, 0.65, 0.8, 0.95};
  const std::vector<mps::Real> exner_half{0.4, 0.575, 0.725, 0.875, 1.0};
  const auto unchanged =
      mps::dry_convective_adjustment(input(stable, mass, exner, exner_half));
  MPS_CHECK(unchanged.adjusted_potential_temperature_k == stable);
  MPS_CHECK_EQ(unchanged.diagnostics.adjusted_layer_count, 0U);

  // Round weights can hide a roundoff-level identity failure, so the same stable column
  // is repeated on the irregular masses and Exner values a hybrid coordinate produces.
  const std::vector<mps::Real> irregular_mass{1080.343469, 186.935824, 576.062394,
                                              1013.056966};
  const std::vector<mps::Real> irregular_exner{0.41, 0.62, 0.81, 0.95};
  const std::vector<mps::Real> irregular_half{0.3, 0.515, 0.715, 0.88, 1.0};
  const auto irregular = mps::dry_convective_adjustment(
      input(stable, irregular_mass, irregular_exner, irregular_half));
  MPS_CHECK(irregular.adjusted_potential_temperature_k == stable);
  MPS_CHECK_EQ(irregular.diagnostics.adjusted_layer_count, 0U);
  MPS_CHECK_EQ(irregular.diagnostics.adjusted_block_count, 0U);
  MPS_CHECK_EQ(irregular.diagnostics.enthalpy_change_j_m2, 0.0);
  MPS_CHECK_EQ(irregular.diagnostics.maximum_temperature_increment_k, 0.0);

  const std::vector<mps::Real> unstable{310.0, 330.0, 290.0, 350.0};
  const auto first =
      mps::dry_convective_adjustment(input(unstable, mass, exner, exner_half));
  const auto second = mps::dry_convective_adjustment(
      input(first.adjusted_potential_temperature_k, mass, exner, exner_half));
  MPS_CHECK(second.adjusted_potential_temperature_k ==
            first.adjusted_potential_temperature_k);
  MPS_CHECK_EQ(second.diagnostics.adjusted_layer_count, 0U);
  MPS_CHECK(first.diagnostics.minimum_theta_difference_after_k >= -1e-10);
  const auto before = enthalpy(unstable, mass, exner);
  const auto after = enthalpy(first.adjusted_potential_temperature_k, mass, exner);
  MPS_CHECK(std::abs(after - before) / before <= 1e-12);
}

MPS_TEST_CASE("dry-energy attribution matches hydrostatic re-diagnosis") {
  const std::vector<mps::Real> theta{310.0, 280.0, 340.0, 300.0};
  const std::vector<mps::Real> mass{1300.0, 2200.0, 3100.0, 3600.0};
  mps::HybridPressureGeometry geometry;
  geometry.exner_full = {0.48, 0.63, 0.78, 0.94};
  geometry.exner_half = {0.40, 0.55, 0.70, 0.86, 1.0};
  const auto result = mps::dry_convective_adjustment(
      input(theta, mass, geometry.exner_full, geometry.exner_half));
  const auto before =
      mps::integrate_hydrostatic_column(geometry, theta, 1004.0, 9.8, 123.0);
  const auto after = mps::integrate_hydrostatic_column(
      geometry, result.adjusted_potential_temperature_k, 1004.0, 9.8, 123.0);
  mps::Real diagnosed_change = 0.0;
  for (std::size_t level = 0; level < theta.size(); ++level) {
    diagnosed_change +=
        mass[level] *
        (717.0 * geometry.exner_full[level] *
             (result.adjusted_potential_temperature_k[level] - theta[level]) +
         after.geopotential_full_m2_s2[level] - before.geopotential_full_m2_s2[level]);
  }
  MPS_CHECK(std::abs(diagnosed_change) > 1e5);
  MPS_CHECK_NEAR(result.diagnostics.dry_energy_attributed_change_j_m2, diagnosed_change,
                 1e-6 * std::abs(diagnosed_change));
}

// A deliberately naive comparator: find the first unstable interface, grow the block
// around it until both neighbours are stable against its weighted mean, flatten it, and
// start over. It is O(K^2) and shares no code with the stack-based merge under test.
[[nodiscard]] std::vector<mps::Real> iterative_reference(
    std::vector<mps::Real> theta, const std::vector<mps::Real>& weight,
    const mps::Real tolerance) {
  for (std::size_t sweep = 0; sweep <= theta.size() * theta.size(); ++sweep) {
    std::size_t seed = theta.size();
    for (std::size_t level = 0; level + 1 < theta.size(); ++level)
      if (theta[level] + tolerance < theta[level + 1]) {
        seed = level;
        break;
      }
    if (seed == theta.size()) return theta;
    std::size_t begin = seed;
    std::size_t end = seed + 2;
    while (true) {
      mps::Real total = 0.0;
      mps::Real weighted = 0.0;
      for (std::size_t level = begin; level < end; ++level) {
        total += weight[level];
        weighted += weight[level] * theta[level];
      }
      const mps::Real mean = weighted / total;
      if (begin > 0 && theta[begin - 1] + tolerance < mean) {
        --begin;
        continue;
      }
      if (end < theta.size() && mean + tolerance < theta[end]) {
        ++end;
        continue;
      }
      for (std::size_t level = begin; level < end; ++level) theta[level] = mean;
      break;
    }
  }
  throw std::runtime_error("iterative reference did not converge");
}

MPS_TEST_CASE("random columns match an independent iterative merger") {
  // A fixed seed and std::mt19937 keep the columns identical across standard libraries.
  std::mt19937 generator(11);
  const auto uniform = [&generator](const mps::Real low, const mps::Real high) {
    return low + (high - low) * static_cast<mps::Real>(generator() - generator.min()) /
                     static_cast<mps::Real>(generator.max() - generator.min());
  };
  std::size_t all_unstable = 0;
  std::size_t thin_layers = 0;
  for (int trial = 0; trial < 200; ++trial) {
    const auto levels = static_cast<std::size_t>(2 + (generator() % 24));
    std::vector<mps::Real> theta(levels), mass(levels), exner_full(levels);
    std::vector<mps::Real> exner_half(levels + 1);
    // Every fourth column is monotonically increasing downward, so the whole column is
    // unstable; every fifth carries a very thin layer against thick neighbours.
    const bool inverted = trial % 4 == 0;
    const bool thin = trial % 5 == 0;
    for (std::size_t level = 0; level < levels; ++level) {
      theta[level] = inverted ? 250.0 + 4.0 * static_cast<mps::Real>(level)
                              : uniform(250.0, 350.0);
      mass[level] =
          (thin && level == levels / 2) ? uniform(0.05, 0.5) : uniform(50.0, 2000.0);
    }
    if (inverted) ++all_unstable;
    if (thin) ++thin_layers;
    exner_half.front() = uniform(0.05, 0.2);
    for (std::size_t level = 0; level < levels; ++level) {
      exner_half[level + 1] = exner_half[level] + uniform(0.01, 0.1);
      exner_full[level] = 0.5 * (exner_half[level] + exner_half[level + 1]);
    }
    std::vector<mps::Real> weight(levels);
    for (std::size_t level = 0; level < levels; ++level)
      weight[level] = 1004.0 * mass[level] * exner_full[level];
    const auto result =
        mps::dry_convective_adjustment(input(theta, mass, exner_full, exner_half));
    const auto expected = iterative_reference(theta, weight, 1e-10);
    for (std::size_t level = 0; level < levels; ++level)
      MPS_CHECK_NEAR(result.adjusted_potential_temperature_k[level], expected[level],
                     1e-9);
    const auto before = enthalpy(theta, mass, exner_full);
    const auto after =
        enthalpy(result.adjusted_potential_temperature_k, mass, exner_full);
    MPS_CHECK(std::abs(after - before) / before <= 1e-12);
    MPS_CHECK(result.diagnostics.minimum_theta_difference_after_k >= -1e-10);
    const auto again = mps::dry_convective_adjustment(
        input(result.adjusted_potential_temperature_k, mass, exner_full, exner_half));
    for (std::size_t level = 0; level < levels; ++level)
      MPS_CHECK(std::abs(again.adjusted_potential_temperature_k[level] -
                         result.adjusted_potential_temperature_k[level]) <= 1e-10);
  }
  MPS_CHECK(all_unstable >= 40);
  MPS_CHECK(thin_layers >= 30);
}

MPS_TEST_CASE("invalid adjustment inputs are rejected") {
  const std::vector<mps::Real> theta{300.0, 310.0};
  const std::vector<mps::Real> mass{1.0, 1.0};
  const std::vector<mps::Real> exner{0.8, 0.95};
  const std::vector<mps::Real> half{0.7, 0.9, 1.0};
  auto invalid = input(theta, mass, exner, half);
  invalid.stability_tolerance_k = -1.0;
  MPS_CHECK_THROWS_AS(mps::dry_convective_adjustment(invalid), std::invalid_argument);
  const std::vector<mps::Real> bad_mass{1.0, 0.0};
  MPS_CHECK_THROWS_AS(
      mps::dry_convective_adjustment(input(theta, bad_mass, exner, half)),
      std::invalid_argument);
}

int main() { return mps::test::run_all(); }
