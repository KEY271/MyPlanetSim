#include "myplanetsim/physics/surface_hydrology.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace mps {
namespace {

void validate_fraction(const Real value, const char* name) {
  if (!std::isfinite(value) || value < 0.0 || value > 1.0)
    throw std::invalid_argument(std::string(name) + " must be in [0, 1]");
}

void validate_bucket(const Real land_fraction, const Real water, const Real capacity) {
  if (!std::isfinite(water) || !std::isfinite(capacity) || capacity < 0.0 ||
      water < 0.0 || water > capacity)
    throw std::invalid_argument("surface bucket state is invalid");
  if (land_fraction == 0.0 && water != 0.0)
    throw std::invalid_argument("all-ocean cells must keep an empty land bucket");
}

}  // namespace

SurfaceEvaporationPartition partition_surface_evaporation(
    const Real saturation_deficit, const Real land_conductance_kg_m2_s,
    const Real ocean_conductance_kg_m2_s, const Real land_fraction,
    const Real initial_land_water_kg_m2, const Real bucket_capacity_kg_m2,
    const Real bucket_wet_threshold_fraction, const Real time_step_s) {
  validate_fraction(land_fraction, "land fraction");
  validate_fraction(bucket_wet_threshold_fraction, "bucket wet threshold fraction");
  validate_bucket(land_fraction, initial_land_water_kg_m2, bucket_capacity_kg_m2);
  if (!std::isfinite(saturation_deficit) || !std::isfinite(land_conductance_kg_m2_s) ||
      !std::isfinite(ocean_conductance_kg_m2_s) || land_conductance_kg_m2_s < 0.0 ||
      ocean_conductance_kg_m2_s < 0.0 || !(time_step_s > 0.0) ||
      !std::isfinite(time_step_s))
    throw std::invalid_argument("surface evaporation arguments are invalid");
  if (land_fraction > 0.0 &&
      (!(bucket_capacity_kg_m2 > 0.0) || !(bucket_wet_threshold_fraction > 0.0)))
    throw std::invalid_argument("land evaporation requires a positive bucket");

  SurfaceEvaporationPartition result;
  result.ocean_flux_kg_m2_s = ocean_conductance_kg_m2_s * saturation_deficit;
  if (land_fraction > 0.0) {
    const Real potential = land_conductance_kg_m2_s * saturation_deficit;
    if (potential < 0.0) {
      result.land_flux_kg_m2_s = potential;
    } else {
      const Real wet_threshold = bucket_wet_threshold_fraction * bucket_capacity_kg_m2;
      if (initial_land_water_kg_m2 - time_step_s * potential >= wet_threshold) {
        result.land_flux_kg_m2_s = potential;
      } else {
        result.land_flux_kg_m2_s = potential * initial_land_water_kg_m2 /
                                   (wet_threshold + time_step_s * potential);
      }
    }
    const Real unconstrained_water =
        initial_land_water_kg_m2 - time_step_s * result.land_flux_kg_m2_s;
    result.runoff_kg_m2_land =
        std::max(0.0, unconstrained_water - bucket_capacity_kg_m2);
    result.final_land_water_kg_m2 =
        std::clamp(unconstrained_water, 0.0, bucket_capacity_kg_m2);
  }
  result.cell_flux_kg_m2_s = land_fraction * result.land_flux_kg_m2_s +
                             (1.0 - land_fraction) * result.ocean_flux_kg_m2_s;
  return result;
}

SurfacePrecipitationUpdate route_surface_precipitation(
    const Real land_fraction, const Real precipitation_kg_m2_cell,
    const Real initial_land_water_kg_m2, const Real bucket_capacity_kg_m2) {
  validate_fraction(land_fraction, "land fraction");
  validate_bucket(land_fraction, initial_land_water_kg_m2, bucket_capacity_kg_m2);
  if (!std::isfinite(precipitation_kg_m2_cell) || precipitation_kg_m2_cell < 0.0)
    throw std::invalid_argument("precipitation must be finite and non-negative");
  if (land_fraction > 0.0 && !(bucket_capacity_kg_m2 > 0.0))
    throw std::invalid_argument("land precipitation requires a positive bucket");

  SurfacePrecipitationUpdate result;
  if (land_fraction == 0.0) {
    result.ocean_gain_kg_m2_cell = precipitation_kg_m2_cell;
    return result;
  }
  const Real unconstrained = initial_land_water_kg_m2 + precipitation_kg_m2_cell;
  result.runoff_kg_m2_land = std::max(0.0, unconstrained - bucket_capacity_kg_m2);
  result.final_land_water_kg_m2 = std::min(unconstrained, bucket_capacity_kg_m2);
  const Real ocean_precipitation = (1.0 - land_fraction) * precipitation_kg_m2_cell;
  const Real cell_runoff = land_fraction * result.runoff_kg_m2_land;
  if (land_fraction < 1.0)
    result.ocean_gain_kg_m2_cell = ocean_precipitation + cell_runoff;
  else
    result.external_outflow_kg_m2_cell = cell_runoff;
  const Real land_gain =
      land_fraction * (result.final_land_water_kg_m2 - initial_land_water_kg_m2);
  result.water_budget_residual_kg_m2 = land_gain + result.ocean_gain_kg_m2_cell +
                                       result.external_outflow_kg_m2_cell -
                                       precipitation_kg_m2_cell;
  return result;
}

}  // namespace mps
