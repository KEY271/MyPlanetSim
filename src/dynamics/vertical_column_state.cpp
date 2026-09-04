#include "myplanetsim/dynamics/vertical_column_state.hpp"

#include <cmath>
#include <stdexcept>

#include "myplanetsim/core/validation.hpp"

namespace mps {

std::vector<Real> flatten_vertical_column_state(const VerticalColumnState& state) {
  const std::size_t nz = state.potential_temperature_mass_k_kg_m2.size();
  if (nz == 0 || state.tracer_mass_kg_m2.size() != nz) {
    throw std::invalid_argument("vertical column state shape is invalid");
  }
  std::vector<Real> flat(1 + 2 * nz);
  flat[0] = state.surface_pressure_pa;
  for (std::size_t k = 0; k < nz; ++k) {
    flat[1 + k] = state.potential_temperature_mass_k_kg_m2[k];
    flat[1 + nz + k] = state.tracer_mass_kg_m2[k];
  }
  return flat;
}

VerticalColumnState unflatten_vertical_column_state(
    const Real time_s, const std::uint64_t step, const std::span<const Real> flat_state,
    const std::size_t levels) {
  if (levels == 0 || flat_state.size() != 1 + 2 * levels) {
    throw std::invalid_argument("vertical column flat state size does not match shape");
  }
  return {.time_s = time_s,
          .step = step,
          .surface_pressure_pa = flat_state[0],
          .potential_temperature_mass_k_kg_m2 = {flat_state.begin() + 1,
                                                 flat_state.begin() + 1 + levels},
          .tracer_mass_kg_m2 = {flat_state.begin() + 1 + levels, flat_state.end()}};
}

void validate_vertical_column_state(const VerticalColumnState& state,
                                    const std::span<const Real> air_mass_kg_m2,
                                    const Real temperature_floor_k,
                                    const std::span<const Real> exner_full,
                                    const bool require_nonnegative_tracer) {
  require_non_negative(state.time_s, "vertical column time");
  require_positive(state.surface_pressure_pa, "vertical column surface pressure");
  require_positive(temperature_floor_k, "vertical temperature floor");
  const std::size_t nz = air_mass_kg_m2.size();
  if (nz == 0 || state.potential_temperature_mass_k_kg_m2.size() != nz ||
      state.tracer_mass_kg_m2.size() != nz || exner_full.size() != nz) {
    throw std::invalid_argument("vertical column state and geometry shapes differ");
  }
  for (std::size_t k = 0; k < nz; ++k) {
    require_positive(air_mass_kg_m2[k], "vertical air mass");
    require_positive(exner_full[k], "vertical full Exner");
    const Real theta = state.potential_temperature_mass_k_kg_m2[k] / air_mass_kg_m2[k];
    if (!std::isfinite(theta) || !(theta * exner_full[k] >= temperature_floor_k)) {
      throw std::runtime_error(
          "vertical column temperature is non-finite or below floor");
    }
    if (!std::isfinite(state.tracer_mass_kg_m2[k]) ||
        (require_nonnegative_tracer && state.tracer_mass_kg_m2[k] < 0.0)) {
      throw std::runtime_error("vertical column tracer mass is invalid");
    }
  }
}

}  // namespace mps
