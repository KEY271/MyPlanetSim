#pragma once

#include "myplanetsim/dynamics/shallow_water_state.hpp"

namespace mps::diagnostics {

struct ShallowWaterInvariants {
  Real mass;
  Real energy;
  Real potential_enstrophy;
  Real axial_angular_momentum;
  Real minimum_depth;
  Real maximum_depth;
  Real minimum_pv;
  Real maximum_pv;
  Real maximum_radial_momentum;
  Real rms_radial_momentum;
};

struct ShallowWaterInvariantRates {
  Real mass;
  Real energy;
  Real potential_enstrophy;
  Real axial_angular_momentum;
};

struct ShallowWaterBudget {
  ShallowWaterInvariantRates flux;
  ShallowWaterInvariantRates coriolis;
  ShallowWaterInvariantRates pressure;
  ShallowWaterInvariantRates diffusion;
  ShallowWaterInvariantRates total;
  ShallowWaterInvariantRates residual;
};

[[nodiscard]] ShallowWaterInvariants diagnose_shallow_water(
    const CubedSphereGrid& grid, const ShallowWaterState& state, Real gravity_m_s2,
    Vec3 rotation_vector_rad_s);
[[nodiscard]] ShallowWaterInvariantRates shallow_water_invariant_rates(
    const CubedSphereGrid& grid, const ShallowWaterState& state,
    const ShallowWaterTendency& tendency, Real gravity_m_s2,
    Vec3 rotation_vector_rad_s);
[[nodiscard]] ShallowWaterBudget make_shallow_water_budget(
    const CubedSphereGrid& grid, const ShallowWaterState& state,
    const ShallowWaterTendency& flux, const ShallowWaterTendency& coriolis,
    const ShallowWaterTendency& pressure, const ShallowWaterTendency& diffusion,
    Real gravity_m_s2, Vec3 rotation_vector_rad_s);

}  // namespace mps::diagnostics
