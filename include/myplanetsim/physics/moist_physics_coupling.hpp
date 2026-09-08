#pragma once

#include "myplanetsim/config/experiment_config.hpp"
#include "myplanetsim/dynamics/dry_hydrostatic_state.hpp"
#include "myplanetsim/dynamics/surface_boundary.hpp"
#include "myplanetsim/grid/cubed_sphere_grid.hpp"
#include "myplanetsim/physics/simple_betts_miller.hpp"

namespace mps {

struct MoistPhysicsStepDiagnostics {
  Real evaporation_kg = 0.0;
  Real convective_precipitation_kg = 0.0;
  Real grid_scale_precipitation_kg = 0.0;
  Real runoff_kg = 0.0;
  Real ocean_water_change_kg = 0.0;
  Real external_outflow_kg = 0.0;
  Real water_budget_residual_kg = 0.0;
  Real moist_enthalpy_budget_residual_j = 0.0;
  Real maximum_relative_humidity = 0.0;
  std::size_t deep_column_count = 0;
  std::size_t shallow_column_count = 0;
  std::size_t inactive_column_count = 0;
};

struct MoistPhysicsCouplingWorkspace {
  SimpleBettsMillerResult convection;
  HybridPressureGeometry vertical_geometry;
};

void apply_moist_column_physics(
    const CubedSphereGrid& grid, const AtmosphericHybridCoordinate& coordinate,
    const SurfaceBoundary& boundary, DryHydrostaticState& state,
    const DryHydrostaticDerived& derived, std::size_t water_vapor_tracer,
    const PlanetParameters& planet, const MoistureParameters& moisture,
    const ConvectionParameters& convection, const SurfaceParameters& surface,
    Real time_step_s, MoistPhysicsStepDiagnostics& diagnostics,
    MoistPhysicsCouplingWorkspace& workspace);

}  // namespace mps
