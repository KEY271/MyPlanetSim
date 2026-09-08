#pragma once

#include <span>
#include <vector>

#include "myplanetsim/dynamics/dry_hydrostatic_state.hpp"
#include "myplanetsim/dynamics/tracer_registry.hpp"
#include "myplanetsim/grid/cubed_sphere_grid.hpp"
#include "myplanetsim/physics/moist_thermodynamics.hpp"

namespace mps {

struct TracerGlobalDiagnostics {
  Real mass_kg = 0.0;
  Real minimum_mixing_ratio = 0.0;
  Real maximum_mixing_ratio = 0.0;
};

struct MoistureGlobalDiagnostics {
  Real atmospheric_water_kg = 0.0;
  Real precipitable_water_kg_m2 = 0.0;
  Real area_mean_relative_humidity = 0.0;
  Real minimum_relative_humidity = 0.0;
  Real maximum_relative_humidity = 0.0;
  Real land_water_kg = 0.0;
  Real area_mean_land_water_kg_m2 = 0.0;
  Real maximum_vapor_mixing_ratio = 0.0;
  Real dilute_limit_exceedance_area_fraction = 0.0;
};

[[nodiscard]] std::vector<TracerGlobalDiagnostics> diagnose_tracers(
    const CubedSphereGrid& grid, const DryHydrostaticState& state,
    const DryHydrostaticDerived& derived, const TracerRegistry& registry);

[[nodiscard]] MoistureGlobalDiagnostics diagnose_moisture(
    const CubedSphereGrid& grid, const DryHydrostaticState& state,
    const DryHydrostaticDerived& derived, std::span<const Real> land_fraction,
    std::size_t water_vapor_tracer, const DiluteMoistThermodynamics& thermodynamics,
    Real dilute_limit_mixing_ratio = 0.1);

}  // namespace mps
