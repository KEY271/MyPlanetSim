#pragma once

#include <span>
#include <vector>

#include "myplanetsim/dynamics/dry_hydrostatic_coupling.hpp"
#include "myplanetsim/dynamics/dry_hydrostatic_fast_modes.hpp"

namespace mps {

// A perturbation in the unchanged dry-hydrostatic prognostic variables. The
// reference momentum is zero, so momentum needs no separate subtraction.
struct DryHydrostaticFastPerturbation {
  std::vector<Real> surface_pressure_pa;
  std::vector<Vec3> horizontal_momentum_mass_kg_m_s;
  std::vector<Real> potential_temperature_mass_k_kg_m2;
};

struct DryHydrostaticFastTendency {
  std::vector<Real> surface_pressure_pa_s;
  DryHydrostaticTransportTendency tendency;
};

// Time-independent coefficients for L_ref. The diagnostic Jacobian is stored
// level-major: d(geopotential[level])/d(M*theta[source_level]).
struct DryHydrostaticFastOperator {
  std::size_t levels = 0;
  std::vector<Real> b_half;
  std::vector<Real> reference_air_mass_kg_m2;
  std::vector<Real> reference_potential_temperature_k;
  std::vector<Real> reference_interface_potential_temperature_k;
  std::vector<Real> reference_specific_volume_m3_kg;
  std::vector<Real> pressure_from_surface_pressure;
  std::vector<Real> geopotential_from_surface_pressure;
  std::vector<Real> geopotential_from_potential_temperature_mass;
};

struct DryHydrostaticFastOperatorWorkspace {
  std::vector<Real> horizontal_air_mass_tendency;
  std::vector<Real> horizontal_potential_temperature_mass_tendency;
  std::vector<Real> pressure_perturbation;
  std::vector<Real> geopotential_perturbation;
  std::vector<Vec3> pressure_gradient;
  std::vector<Vec3> geopotential_gradient;
  VerticalMassFlux vertical_mass_flux;
};

[[nodiscard]] DryHydrostaticFastOperator make_dry_hydrostatic_fast_operator(
    const AtmosphericHybridCoordinate& coordinate, const PlanetParameters& planet,
    const DryHydrostaticReferenceColumn& reference);

[[nodiscard]] DryHydrostaticFastPerturbation make_dry_hydrostatic_fast_perturbation(
    const DryHydrostaticState& state, const DryHydrostaticReferenceColumn& reference);

void apply_dry_hydrostatic_fast_operator(
    const CubedSphereGrid& grid, const PlanetParameters& planet,
    const DryHydrostaticFastOperator& fast_operator,
    const DryHydrostaticFastPerturbation& perturbation,
    DryHydrostaticFastTendency& result, DryHydrostaticFastOperatorWorkspace& workspace);

[[nodiscard]] DryHydrostaticFastTendency subtract_dry_hydrostatic_fast_tendency(
    const DryHydrostaticFastTendency& left, const DryHydrostaticFastTendency& right);
[[nodiscard]] DryHydrostaticFastTendency add_dry_hydrostatic_fast_tendency(
    const DryHydrostaticFastTendency& left, const DryHydrostaticFastTendency& right);

}  // namespace mps
