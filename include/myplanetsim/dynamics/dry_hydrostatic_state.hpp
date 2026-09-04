#pragma once

#include <cstdint>
#include <span>
#include <string_view>
#include <vector>

#include "myplanetsim/core/planet_parameters.hpp"
#include "myplanetsim/geometry/vec3.hpp"
#include "myplanetsim/vertical/hybrid_pressure_coordinate.hpp"
#include "myplanetsim/vertical/hydrostatic_column.hpp"

namespace mps {

inline constexpr std::string_view kDryHydrostaticCheckpointLayout =
    "dry_hydrostatic_cell_column_v1";

struct DryHydrostaticState {
  Real time_s = 0.0;
  std::uint64_t step = 0;
  std::vector<Real> surface_pressure_pa;
  std::vector<Vec3> horizontal_momentum_mass_kg_m_s;
  std::vector<Real> potential_temperature_mass_k_kg_m2;
  std::vector<Real> tracer_mass_kg_m2;
};

struct DryHydrostaticDerived {
  std::size_t cells = 0;
  std::size_t levels = 0;
  std::vector<Real> pressure_pa;
  std::vector<Real> air_mass_kg_m2;
  std::vector<Vec3> velocity_m_s;
  std::vector<Real> potential_temperature_k;
  std::vector<Real> tracer_mixing_ratio;
  std::vector<Real> temperature_k;
  std::vector<Real> geopotential_m2_s2;
};

[[nodiscard]] constexpr std::size_t dry_hydrostatic_offset(
    const std::size_t cell, const std::size_t level,
    const std::size_t levels) noexcept {
  return cell * levels + level;
}
[[nodiscard]] std::vector<Real> flatten_dry_hydrostatic_state(
    const DryHydrostaticState& state, std::size_t levels);
[[nodiscard]] DryHydrostaticState unflatten_dry_hydrostatic_state(
    Real time_s, std::uint64_t step, std::span<const Real> values, std::size_t cells,
    std::size_t levels);
[[nodiscard]] DryHydrostaticDerived diagnose_dry_hydrostatic_state(
    const DryHydrostaticState& state, const AtmosphericHybridCoordinate& coordinate,
    const PlanetParameters& planet);
void validate_dry_hydrostatic_state(const DryHydrostaticState& state,
                                    const DryHydrostaticDerived& derived,
                                    std::span<const Vec3> cell_centres,
                                    Real minimum_surface_pressure_pa,
                                    Real maximum_surface_pressure_pa,
                                    Real temperature_floor_k,
                                    bool require_nonnegative_tracer = true);

}  // namespace mps
