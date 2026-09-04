#pragma once

#include <cstdint>
#include <span>
#include <string_view>
#include <vector>

#include "myplanetsim/core/types.hpp"

namespace mps {

inline constexpr std::string_view kVerticalColumnCheckpointLayout =
    "vertical_column_hybrid_v1";

struct VerticalColumnState {
  Real time_s = 0.0;
  std::uint64_t step = 0;
  Real surface_pressure_pa = 0.0;
  std::vector<Real> potential_temperature_mass_k_kg_m2;
  std::vector<Real> tracer_mass_kg_m2;
};

[[nodiscard]] std::vector<Real> flatten_vertical_column_state(
    const VerticalColumnState& state);
[[nodiscard]] VerticalColumnState unflatten_vertical_column_state(
    Real time_s, std::uint64_t step, std::span<const Real> flat_state,
    std::size_t levels);
void validate_vertical_column_state(const VerticalColumnState& state,
                                    std::span<const Real> air_mass_kg_m2,
                                    Real temperature_floor_k,
                                    std::span<const Real> exner_full,
                                    bool require_nonnegative_tracer = true);

}  // namespace mps
