#include "myplanetsim/physics/dry_convective_adjustment.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "myplanetsim/core/validation.hpp"

namespace mps {
namespace {

struct StabilitySummary {
  Real minimum_difference_k = 0.0;
  Real unstable_fraction = 0.0;
};

[[nodiscard]] StabilitySummary summarize_stability(const std::span<const Real> theta,
                                                   const Real tolerance_k) {
  if (theta.size() < 2) return {};
  Real minimum = std::numeric_limits<Real>::infinity();
  std::size_t unstable = 0;
  for (std::size_t level = 0; level + 1 < theta.size(); ++level) {
    const Real difference = theta[level] - theta[level + 1];
    minimum = std::min(minimum, difference);
    if (difference < -tolerance_k) ++unstable;
  }
  return {.minimum_difference_k = minimum,
          .unstable_fraction =
              static_cast<Real>(unstable) / static_cast<Real>(theta.size() - 1)};
}

}  // namespace

void dry_convective_adjustment(const DryConvectiveAdjustmentInput& input,
                               DryConvectiveAdjustmentResult& result,
                               DryConvectiveAdjustmentWorkspace& workspace) {
  const std::size_t levels = input.potential_temperature_k.size();
  if (levels == 0 || input.air_mass_kg_m2.size() != levels ||
      input.exner_full.size() != levels || input.exner_half.size() != levels + 1)
    throw std::invalid_argument("dry convective adjustment shape mismatch");
  require_positive(input.heat_capacity_cp_j_kg_k, "heat capacity cp");
  require_positive(input.heat_capacity_cv_j_kg_k, "heat capacity cv");
  require_non_negative(input.stability_tolerance_k, "stability tolerance");

  result.adjusted_potential_temperature_k.resize(levels);
  result.potential_temperature_mass_increment_k_kg_m2.resize(levels);
  workspace.blocks.clear();
  workspace.blocks.reserve(levels);

  for (std::size_t level = 0; level < levels; ++level) {
    require_positive(input.potential_temperature_k[level], "potential temperature");
    require_positive(input.air_mass_kg_m2[level], "air mass");
    require_positive(input.exner_full[level], "full-level Exner function");
    require_positive(input.exner_half[level], "half-level Exner function");
    const Real weight = input.heat_capacity_cp_j_kg_k * input.air_mass_kg_m2[level] *
                        input.exner_full[level];
    require_positive(weight, "convective adjustment weight");
    workspace.blocks.push_back(
        {.begin = level,
         .end = level + 1,
         .weight = weight,
         .weighted_theta = weight * input.potential_temperature_k[level]});
    while (workspace.blocks.size() >= 2) {
      const auto& upper = workspace.blocks[workspace.blocks.size() - 2];
      const auto& lower = workspace.blocks.back();
      const Real upper_mean = upper.weighted_theta / upper.weight;
      const Real lower_mean = lower.weighted_theta / lower.weight;
      if (upper_mean + input.stability_tolerance_k >= lower_mean) break;
      auto merged = DryConvectiveAdjustmentWorkspace::Block{
          .begin = upper.begin,
          .end = lower.end,
          .weight = upper.weight + lower.weight,
          .weighted_theta = upper.weighted_theta + lower.weighted_theta,
      };
      workspace.blocks.pop_back();
      workspace.blocks.back() = merged;
    }
  }
  require_positive(input.exner_half.back(), "half-level Exner function");

  std::size_t adjusted_blocks = 0;
  for (const auto& block : workspace.blocks) {
    if (block.end - block.begin == 1) {
      // A layer that never merged is copied rather than recovered from its own weighted
      // mean: dividing the weight back out perturbs it at roundoff, which would make a
      // stable column report a spurious adjustment and break idempotence.
      result.adjusted_potential_temperature_k[block.begin] =
          input.potential_temperature_k[block.begin];
      continue;
    }
    ++adjusted_blocks;
    const Real mean = block.weighted_theta / block.weight;
    for (std::size_t level = block.begin; level < block.end; ++level)
      result.adjusted_potential_temperature_k[level] = mean;
  }

  const auto before =
      summarize_stability(input.potential_temperature_k, input.stability_tolerance_k);
  const auto after = summarize_stability(result.adjusted_potential_temperature_k,
                                         input.stability_tolerance_k);
  Real enthalpy_change = 0.0;
  Real maximum_temperature_increment = 0.0;
  std::size_t adjusted_layers = 0;
  for (std::size_t level = 0; level < levels; ++level) {
    const Real delta_theta = result.adjusted_potential_temperature_k[level] -
                             input.potential_temperature_k[level];
    result.potential_temperature_mass_increment_k_kg_m2[level] =
        input.air_mass_kg_m2[level] * delta_theta;
    enthalpy_change += input.heat_capacity_cp_j_kg_k * input.exner_full[level] *
                       result.potential_temperature_mass_increment_k_kg_m2[level];
    maximum_temperature_increment = std::max(
        maximum_temperature_increment, std::abs(input.exner_full[level] * delta_theta));
    if (delta_theta != 0.0) ++adjusted_layers;
  }

  Real dry_energy_change = 0.0;
  Real delta_geopotential_half = 0.0;
  for (std::size_t reverse = levels; reverse > 0; --reverse) {
    const std::size_t level = reverse - 1;
    const Real delta_theta = result.adjusted_potential_temperature_k[level] -
                             input.potential_temperature_k[level];
    const Real delta_geopotential_full =
        delta_geopotential_half +
        input.heat_capacity_cp_j_kg_k * delta_theta *
            (input.exner_half[level + 1] - input.exner_full[level]);
    dry_energy_change +=
        input.air_mass_kg_m2[level] *
        (input.heat_capacity_cv_j_kg_k * input.exner_full[level] * delta_theta +
         delta_geopotential_full);
    delta_geopotential_half += input.heat_capacity_cp_j_kg_k * delta_theta *
                               (input.exner_half[level + 1] - input.exner_half[level]);
  }

  result.diagnostics = {
      .minimum_theta_difference_before_k = before.minimum_difference_k,
      .minimum_theta_difference_after_k = after.minimum_difference_k,
      .unstable_interface_fraction_before = before.unstable_fraction,
      .unstable_interface_fraction_after = after.unstable_fraction,
      .adjusted_layer_count = adjusted_layers,
      .adjusted_block_count = adjusted_blocks,
      .maximum_temperature_increment_k = maximum_temperature_increment,
      .enthalpy_change_j_m2 = enthalpy_change,
      .dry_energy_attributed_change_j_m2 = dry_energy_change,
  };
}

DryConvectiveAdjustmentResult dry_convective_adjustment(
    const DryConvectiveAdjustmentInput& input) {
  DryConvectiveAdjustmentResult result;
  DryConvectiveAdjustmentWorkspace workspace;
  dry_convective_adjustment(input, result, workspace);
  return result;
}

}  // namespace mps
