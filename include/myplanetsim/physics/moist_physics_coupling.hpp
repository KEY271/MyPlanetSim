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
  Real maximum_temperature_increment_k = 0.0;
  Real maximum_vapor_increment = 0.0;
  Real cape_area_time_integral_j_m2_s_kg = 0.0;
  Real cin_area_time_integral_j_m2_s_kg = 0.0;
  Real lcl_pressure_area_time_integral_pa_m2_s = 0.0;
  Real convection_top_pressure_area_time_integral_pa_m2_s = 0.0;
  Real active_area_time_m2_s = 0.0;
  Real deep_area_time_m2_s = 0.0;
  Real shallow_area_time_m2_s = 0.0;
  Real inactive_area_time_m2_s = 0.0;
  Real model_top_area_time_m2_s = 0.0;
  std::size_t deep_column_count = 0;
  std::size_t shallow_column_count = 0;
  std::size_t inactive_column_count = 0;
};

struct MoistPhysicsCouplingWorkspace {
  SimpleBettsMillerReference convection_reference;
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
