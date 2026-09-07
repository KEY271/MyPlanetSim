#pragma once

#include "myplanetsim/core/types.hpp"

namespace mps {

struct SurfaceEvaporationPartition {
  Real land_flux_kg_m2_s = 0.0;
  Real ocean_flux_kg_m2_s = 0.0;
  Real cell_flux_kg_m2_s = 0.0;
  Real final_land_water_kg_m2 = 0.0;
  Real runoff_kg_m2_land = 0.0;
};

[[nodiscard]] SurfaceEvaporationPartition partition_surface_evaporation(
    Real saturation_deficit, Real land_conductance_kg_m2_s,
    Real ocean_conductance_kg_m2_s, Real land_fraction, Real initial_land_water_kg_m2,
    Real bucket_capacity_kg_m2, Real bucket_wet_threshold_fraction, Real time_step_s);

struct SurfacePrecipitationUpdate {
  Real final_land_water_kg_m2 = 0.0;
  Real runoff_kg_m2_land = 0.0;
  Real ocean_gain_kg_m2_cell = 0.0;
  Real external_outflow_kg_m2_cell = 0.0;
  Real water_budget_residual_kg_m2 = 0.0;
};

[[nodiscard]] SurfacePrecipitationUpdate route_surface_precipitation(
    Real land_fraction, Real precipitation_kg_m2_cell, Real initial_land_water_kg_m2,
    Real bucket_capacity_kg_m2);

}  // namespace mps
