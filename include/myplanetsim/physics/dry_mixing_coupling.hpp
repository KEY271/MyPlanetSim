#pragma once

#include "myplanetsim/config/experiment_config.hpp"
#include "myplanetsim/dynamics/dry_hydrostatic_state.hpp"
#include "myplanetsim/dynamics/surface_boundary.hpp"
#include "myplanetsim/grid/cubed_sphere_grid.hpp"
#include "myplanetsim/physics/boundary_layer.hpp"

namespace mps {

struct DryMixingStepDiagnostics {
  Real sensible_to_atmosphere_energy_j = 0.0;
  Real surface_stress_impulse_magnitude_n_s = 0.0;
  Real mean_boundary_layer_height_m = 0.0;
  Real maximum_momentum_diffusivity_m2_s = 0.0;
  Real maximum_heat_diffusivity_m2_s = 0.0;
  Real shallow_unresolved_area_fraction = 0.0;
  Real model_top_area_fraction = 0.0;
  Real atmospheric_heat_change_j = 0.0;
  Real surface_heat_change_j = 0.0;
  Real physical_shear_dissipation_j = 0.0;
  Real physical_surface_drag_dissipation_j = 0.0;
  Real backward_euler_dissipation_j = 0.0;
  Real returned_dissipation_heat_j = 0.0;
  Real heat_budget_residual_j = 0.0;
  Real kinetic_energy_change_j = 0.0;
  Real kinetic_energy_identity_residual_j = 0.0;
  Real momentum_budget_residual_n_s = 0.0;
  Real tracer_mass_change_kg = 0.0;
};

struct DryMixingCouplingWorkspace {
  HybridPressureGeometry vertical_geometry;
  HydrostaticColumn hydrostatic_column;
  std::vector<Real> height_half_m;
  std::vector<Real> height_full_m;
  BoundaryLayerBulkResult bulk;
  BoundaryLayerColumnResult column;
  BoundaryLayerColumnWorkspace column_workspace;
};

void apply_dry_boundary_layer(
    const CubedSphereGrid& grid, const AtmosphericHybridCoordinate& coordinate,
    const SurfaceBoundary& boundary, DryHydrostaticState& state,
    const DryHydrostaticDerived& derived, const PlanetParameters& planet,
    const SurfaceParameters& surface, const BoundaryLayerParameters& boundary_layer,
    Real time_step_s, DryMixingStepDiagnostics& diagnostics,
    DryMixingCouplingWorkspace& workspace);

}  // namespace mps
