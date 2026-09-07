#pragma once

#include <cstddef>
#include <span>
#include <vector>

#include "myplanetsim/core/types.hpp"

namespace mps {

struct DryConvectiveAdjustmentInput {
  std::span<const Real> potential_temperature_k;
  std::span<const Real> air_mass_kg_m2;
  std::span<const Real> exner_full;
  std::span<const Real> exner_half;
  Real heat_capacity_cp_j_kg_k = 0.0;
  Real heat_capacity_cv_j_kg_k = 0.0;
  Real stability_tolerance_k = 1e-10;
};

struct DryConvectiveAdjustmentDiagnostics {
  Real minimum_theta_difference_before_k = 0.0;
  Real minimum_theta_difference_after_k = 0.0;
  Real unstable_interface_fraction_before = 0.0;
  Real unstable_interface_fraction_after = 0.0;
  std::size_t adjusted_layer_count = 0;
  std::size_t adjusted_block_count = 0;
  Real maximum_temperature_increment_k = 0.0;
  Real enthalpy_change_j_m2 = 0.0;
  Real dry_energy_attributed_change_j_m2 = 0.0;
};

struct DryConvectiveAdjustmentResult {
  std::vector<Real> adjusted_potential_temperature_k;
  std::vector<Real> potential_temperature_mass_increment_k_kg_m2;
  DryConvectiveAdjustmentDiagnostics diagnostics;
};

struct DryConvectiveAdjustmentWorkspace {
  struct Block {
    std::size_t begin = 0;
    std::size_t end = 0;
    Real weight = 0.0;
    Real weighted_theta = 0.0;
  };
  std::vector<Block> blocks;
};

void dry_convective_adjustment(const DryConvectiveAdjustmentInput& input,
                               DryConvectiveAdjustmentResult& result,
                               DryConvectiveAdjustmentWorkspace& workspace);

[[nodiscard]] DryConvectiveAdjustmentResult dry_convective_adjustment(
    const DryConvectiveAdjustmentInput& input);

}  // namespace mps
