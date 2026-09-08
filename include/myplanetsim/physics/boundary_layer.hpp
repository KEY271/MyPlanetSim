#pragma once

#include <span>
#include <vector>

#include "myplanetsim/core/types.hpp"
#include "myplanetsim/geometry/vec3.hpp"
#include "myplanetsim/physics/moist_thermodynamics.hpp"

namespace mps {

struct SurfaceRoughness {
  Real momentum_m = 0.0;
  Real heat_m = 0.0;
};

struct BoundaryLayerBulkInput {
  std::span<const Real> potential_temperature_k;
  std::span<const Real> temperature_k;
  std::span<const Vec3> velocity_m_s;
  std::span<const Real> pressure_half_pa;
  std::span<const Real> exner_half;
  std::span<const Real> height_half_m;
  std::span<const Real> height_full_m;
  Real surface_temperature_k = 0.0;
  Real surface_exner = 0.0;
  Real gravity_m_s2 = 0.0;
  Real gas_constant_j_kg_k = 0.0;
  Real heat_capacity_cp_j_kg_k = 0.0;
  Real critical_richardson = 1.0;
  Real turbulent_prandtl = 1.0;
  Real gustiness_m_s = 1.0;
  Real land_fraction = 0.0;
  SurfaceRoughness land_roughness;
  SurfaceRoughness ocean_roughness;
};

struct BoundaryLayerBulkDiagnostics {
  Real surface_richardson = 0.0;
  Real surface_stability_factor = 0.0;
  Real drag_coefficient = 0.0;
  Real heat_exchange_coefficient = 0.0;
  Real boundary_layer_height_m = 0.0;
  Real maximum_momentum_diffusivity_m2_s = 0.0;
  Real maximum_heat_diffusivity_m2_s = 0.0;
  Real sensible_heat_flux_w_m2 = 0.0;
  Vec3 surface_stress_kg_m_s2{};
  bool shallow_stable_layer_unresolved = false;
  bool reaches_model_top = false;
};

struct BoundaryLayerBulkResult {
  std::vector<Real> density_half_kg_m3;
  std::vector<Real> eddy_diffusivity_momentum_m2_s;
  std::vector<Real> eddy_diffusivity_heat_m2_s;
  std::vector<Real> eddy_diffusivity_tracer_m2_s;
  Real surface_heat_conductance_w_m2_k = 0.0;
  Real surface_drag_conductance_kg_m2_s = 0.0;
  Real surface_water_conductance_land_kg_m2_s = 0.0;
  Real surface_water_conductance_ocean_kg_m2_s = 0.0;
  BoundaryLayerBulkDiagnostics diagnostics;
};

void diagnose_boundary_layer_column(const BoundaryLayerBulkInput& input,
                                    BoundaryLayerBulkResult& result);

[[nodiscard]] BoundaryLayerBulkResult diagnose_boundary_layer_column(
    const BoundaryLayerBulkInput& input);

struct BoundaryLayerColumnInput {
  std::span<const Real> potential_temperature_k;
  std::span<const Vec3> velocity_m_s;
  std::span<const Real> tracer_mixing_ratio;
  std::span<const Real> air_mass_kg_m2;
  std::span<const Real> exner_full;
  std::span<const Real> exner_half;
  std::span<const Real> height_full_m;
  std::span<const Real> density_half_kg_m3;
  std::span<const Real> eddy_diffusivity_momentum_m2_s;
  std::span<const Real> eddy_diffusivity_heat_m2_s;
  std::span<const Real> eddy_diffusivity_tracer_m2_s;
  Real surface_temperature_k = 0.0;
  Real surface_exner = 0.0;
  Real surface_heat_capacity_j_m2_k = 0.0;
  Real surface_heat_conductance_w_m2_k = 0.0;
  Real surface_drag_conductance_kg_m2_s = 0.0;
  Real heat_capacity_cp_j_kg_k = 0.0;
  Real time_step_s = 0.0;
  bool enable_surface_water_exchange = false;
  Real surface_pressure_pa = 0.0;
  Real land_fraction = 0.0;
  Real land_water_kg_m2 = 0.0;
  Real bucket_capacity_kg_m2 = 0.0;
  Real bucket_wet_threshold_fraction = 0.75;
  Real surface_water_conductance_land_kg_m2_s = 0.0;
  Real surface_water_conductance_ocean_kg_m2_s = 0.0;
  DiluteMoistThermodynamics moist_thermodynamics{};
};

struct BoundaryLayerColumnDiagnostics {
  Real atmospheric_heat_change_j_m2 = 0.0;
  Real surface_heat_change_j_m2 = 0.0;
  Real physical_shear_dissipation_j_m2 = 0.0;
  Real physical_surface_drag_dissipation_j_m2 = 0.0;
  Real backward_euler_dissipation_j_m2 = 0.0;
  Real returned_dissipation_heat_j_m2 = 0.0;
  Real heat_budget_residual_j_m2 = 0.0;
  Real kinetic_energy_change_j_m2 = 0.0;
  Real kinetic_energy_identity_residual_j_m2 = 0.0;
  Vec3 atmospheric_momentum_change_kg_m_s{};
  Vec3 surface_stress_impulse_kg_m_s{};
  Vec3 momentum_budget_residual_kg_m_s{};
  Real tracer_mass_change_kg_m2 = 0.0;
  Real evaporation_kg_m2 = 0.0;
  Real land_evaporation_kg_m2_land = 0.0;
  Real ocean_evaporation_kg_m2_ocean = 0.0;
  Real runoff_kg_m2_land = 0.0;
  Real ocean_water_change_kg_m2 = 0.0;
  Real external_outflow_kg_m2 = 0.0;
  Real latent_surface_heat_change_j_m2 = 0.0;
  Real moist_enthalpy_budget_residual_j_m2 = 0.0;
  Real water_budget_residual_kg_m2 = 0.0;
  std::size_t surface_water_iterations = 0;
};

struct BoundaryLayerColumnResult {
  std::vector<Real> potential_temperature_k;
  std::vector<Vec3> velocity_m_s;
  std::vector<Real> tracer_mixing_ratio;
  Real surface_temperature_k = 0.0;
  std::vector<Real> heat_flux_w_m2;
  std::vector<Vec3> momentum_flux_kg_m_s2;
  std::vector<Real> tracer_flux_kg_m2_s;
  std::vector<Real> dissipated_heat_j_m2;
  Real land_water_kg_m2 = 0.0;
  Real surface_water_flux_kg_m2_s = 0.0;
  Real land_water_flux_kg_m2_s = 0.0;
  Real ocean_water_flux_kg_m2_s = 0.0;
  BoundaryLayerColumnDiagnostics diagnostics;
};

struct BoundaryLayerColumnWorkspace {
  std::vector<Real> lower;
  std::vector<Real> diagonal;
  std::vector<Real> upper;
  std::vector<Real> rhs;
  std::vector<Real> solution;
  std::vector<Real> momentum_conductance;
  std::vector<Real> heat_conductance;
  std::vector<Real> tracer_conductance;
  std::vector<Real> closed_momentum_conductance;
  std::vector<Real> momentum_capacity;
  std::vector<Real> heat_capacity;
  std::vector<Real> tracer_capacity;
  std::vector<Real> heat_increment;
  std::vector<Real> tracer_surface_response;
  std::vector<Real> heat_surface_response;
};

void implicit_boundary_layer_column(const BoundaryLayerColumnInput& input,
                                    BoundaryLayerColumnResult& result,
                                    BoundaryLayerColumnWorkspace& workspace);

[[nodiscard]] BoundaryLayerColumnResult implicit_boundary_layer_column(
    const BoundaryLayerColumnInput& input);

}  // namespace mps
