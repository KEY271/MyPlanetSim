#pragma once

#include <span>
#include <vector>

#include "myplanetsim/core/planet_parameters.hpp"
#include "myplanetsim/dynamics/dry_hydrostatic_state.hpp"
#include "myplanetsim/grid/cubed_sphere_grid.hpp"

namespace mps {

struct HeldSuarezRates {
  Real equilibrium_temperature_k;
  Real temperature_relaxation_rate_s_1;
  Real rayleigh_drag_rate_s_1;
};

struct HeldSuarezDiagnostics {
  Real potential_temperature_mass_rate_k_kg_s;
  Real eastward_momentum_rate_n;
  Real northward_momentum_rate_n;
  Real thermal_energy_rate_w;
  Real rayleigh_drag_work_w;
};

struct HeldSuarezTendency {
  std::vector<Real> potential_temperature_mass_k_kg_m2_s;
  std::vector<Vec3> horizontal_momentum_mass_kg_m_s2;
  std::vector<Real> equilibrium_temperature_k;
  std::vector<Real> temperature_relaxation_rate_s_1;
  std::vector<Real> rayleigh_drag_rate_s_1;
  HeldSuarezDiagnostics diagnostics{};
};
struct HeldSuarezWorkspace {
  HybridPressureGeometry vertical_geometry;
  std::vector<Real> theta_mass_rates;
};

[[nodiscard]] HeldSuarezRates held_suarez_rates(Real latitude_rad, Real pressure_pa,
                                                Real surface_pressure_pa,
                                                const PlanetParameters& planet);

[[nodiscard]] HeldSuarezTendency held_suarez_tendency(
    const CubedSphereGrid& grid, const AtmosphericHybridCoordinate& coordinate,
    const DryHydrostaticDerived& derived, std::span<const Real> surface_pressure_pa,
    const PlanetParameters& planet);
void held_suarez_tendency(const CubedSphereGrid& grid,
                          const AtmosphericHybridCoordinate& coordinate,
                          const DryHydrostaticDerived& derived,
                          std::span<const Real> surface_pressure_pa,
                          const PlanetParameters& planet, HeldSuarezTendency& result,
                          HeldSuarezWorkspace& workspace);

}  // namespace mps
