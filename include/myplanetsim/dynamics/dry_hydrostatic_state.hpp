#pragma once

#include <cstdint>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#include "myplanetsim/core/planet_parameters.hpp"
#include "myplanetsim/dynamics/tracer_registry.hpp"
#include "myplanetsim/geometry/vec3.hpp"
#include "myplanetsim/vertical/hybrid_pressure_coordinate.hpp"
#include "myplanetsim/vertical/hydrostatic_column.hpp"

namespace mps {

inline constexpr std::string_view kDryHydrostaticCheckpointLayout =
    "dry_hydrostatic_cell_column_v1";
inline constexpr std::string_view kDryHydrostaticSurfaceCheckpointLayout =
    "dry_hydrostatic_surface_v1";
inline constexpr std::string_view kDryHydrostaticMultitracerCheckpointLayout =
    "dry_hydrostatic_multitracer_v1";
inline constexpr std::string_view kDryHydrostaticMultitracerSurfaceCheckpointLayout =
    "dry_hydrostatic_multitracer_surface_v1";
inline constexpr std::string_view kDryHydrostaticMoistCheckpointLayout =
    "dry_hydrostatic_moist_v1";

struct DryHydrostaticState {
  Real time_s = 0.0;
  std::uint64_t step = 0;
  std::size_t tracer_count = 1;
  std::vector<Real> surface_pressure_pa;
  std::vector<Vec3> horizontal_momentum_mass_kg_m_s;
  std::vector<Real> potential_temperature_mass_k_kg_m2;
  std::vector<Real> tracer_mass_kg_m2;
  std::vector<Real> surface_temperature_k;
  std::vector<Real> land_water_kg_m2;
  Real cumulative_convective_precipitation_kg = 0.0;
  Real cumulative_grid_scale_precipitation_kg = 0.0;
  Real cumulative_evaporation_kg = 0.0;
  Real cumulative_runoff_kg = 0.0;
  Real cumulative_ocean_water_change_kg = 0.0;
  Real cumulative_external_outflow_kg = 0.0;
};

struct DryHydrostaticDerived {
  std::size_t cells = 0;
  std::size_t levels = 0;
  std::size_t tracer_count = 1;
  std::vector<Real> pressure_pa;
  std::vector<Real> exner_half;
  std::vector<Real> exner_full;
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
[[nodiscard]] constexpr std::size_t dry_hydrostatic_tracer_offset(
    const std::size_t tracer, const std::size_t cell, const std::size_t level,
    const std::size_t cells, const std::size_t levels) noexcept {
  return (tracer * cells + cell) * levels + level;
}
[[nodiscard]] std::span<Real> dry_hydrostatic_tracer_component(
    DryHydrostaticState& state, std::size_t tracer, std::size_t cells,
    std::size_t levels);
[[nodiscard]] std::span<const Real> dry_hydrostatic_tracer_component(
    const DryHydrostaticState& state, std::size_t tracer, std::size_t cells,
    std::size_t levels);
[[nodiscard]] std::vector<Real> flatten_dry_hydrostatic_state(
    const DryHydrostaticState& state, std::size_t levels);
[[nodiscard]] DryHydrostaticState unflatten_dry_hydrostatic_state(
    Real time_s, std::uint64_t step, std::span<const Real> values, std::size_t cells,
    std::size_t levels);
[[nodiscard]] std::vector<Real> flatten_dry_hydrostatic_surface_state(
    const DryHydrostaticState& state, std::size_t levels);
[[nodiscard]] DryHydrostaticState unflatten_dry_hydrostatic_surface_state(
    Real time_s, std::uint64_t step, std::span<const Real> values, std::size_t cells,
    std::size_t levels);
[[nodiscard]] std::string dry_hydrostatic_multitracer_checkpoint_layout(
    const TracerRegistry& registry, bool include_surface);
[[nodiscard]] std::vector<Real> flatten_dry_hydrostatic_multitracer_state(
    const DryHydrostaticState& state, std::size_t levels, bool include_surface);
[[nodiscard]] DryHydrostaticState unflatten_dry_hydrostatic_multitracer_state(
    Real time_s, std::uint64_t step, std::span<const Real> values, std::size_t cells,
    std::size_t levels, std::size_t tracer_count, bool include_surface);
[[nodiscard]] std::string dry_hydrostatic_moist_checkpoint_layout(
    const TracerRegistry& registry);
[[nodiscard]] std::vector<Real> flatten_dry_hydrostatic_moist_state(
    const DryHydrostaticState& state, std::size_t levels);
[[nodiscard]] DryHydrostaticState unflatten_dry_hydrostatic_moist_state(
    Real time_s, std::uint64_t step, std::span<const Real> values, std::size_t cells,
    std::size_t levels, std::size_t tracer_count);
[[nodiscard]] DryHydrostaticDerived diagnose_dry_hydrostatic_state(
    const DryHydrostaticState& state, const AtmosphericHybridCoordinate& coordinate,
    const PlanetParameters& planet);
[[nodiscard]] DryHydrostaticDerived diagnose_dry_hydrostatic_state(
    const DryHydrostaticState& state, const AtmosphericHybridCoordinate& coordinate,
    const PlanetParameters& planet, std::span<const Real> surface_geopotential_m2_s2);
void diagnose_dry_hydrostatic_state(
    const DryHydrostaticState& state, const AtmosphericHybridCoordinate& coordinate,
    const PlanetParameters& planet, std::span<const Real> surface_geopotential_m2_s2,
    DryHydrostaticDerived& result, HybridPressureGeometry& geometry_workspace,
    HydrostaticColumn& column_workspace, std::span<Real> theta_workspace);
void validate_dry_hydrostatic_state(const DryHydrostaticState& state,
                                    const DryHydrostaticDerived& derived,
                                    std::span<const Vec3> cell_centres,
                                    Real minimum_surface_pressure_pa,
                                    Real maximum_surface_pressure_pa,
                                    Real temperature_floor_k,
                                    bool require_nonnegative_tracer = true);

}  // namespace mps
