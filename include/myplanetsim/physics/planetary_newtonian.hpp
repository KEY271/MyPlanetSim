#pragma once

#include "myplanetsim/config/experiment_config.hpp"
#include "myplanetsim/core/orbit_parameters.hpp"
#include "myplanetsim/physics/held_suarez.hpp"

namespace mps {

[[nodiscard]] HeldSuarezRates planetary_newtonian_rates(Real latitude_rad, Real chi,
                                                        Real pressure_pa,
                                                        Real surface_pressure_pa,
                                                        const PlanetParameters& planet);

[[nodiscard]] HeldSuarezTendency planetary_newtonian_tendency(
    const CubedSphereGrid& grid, const AtmosphericHybridCoordinate& coordinate,
    const DryHydrostaticDerived& derived, std::span<const Real> surface_pressure_pa,
    const PlanetParameters& planet, ForcingGeometry geometry,
    const OrbitState* orbit_state = nullptr);
void planetary_newtonian_tendency(
    const CubedSphereGrid& grid, const AtmosphericHybridCoordinate& coordinate,
    const DryHydrostaticDerived& derived, std::span<const Real> surface_pressure_pa,
    const PlanetParameters& planet, ForcingGeometry geometry,
    const OrbitState* orbit_state, HeldSuarezTendency& result,
    HeldSuarezWorkspace& workspace);

}  // namespace mps
