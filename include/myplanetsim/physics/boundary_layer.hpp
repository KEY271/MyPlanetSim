#pragma once

#include <span>
#include <vector>

#include "myplanetsim/core/types.hpp"
#include "myplanetsim/geometry/vec3.hpp"

namespace mps {

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
};

void implicit_boundary_layer_column(const BoundaryLayerColumnInput& input,
                                    BoundaryLayerColumnResult& result,
                                    BoundaryLayerColumnWorkspace& workspace);

[[nodiscard]] BoundaryLayerColumnResult implicit_boundary_layer_column(
    const BoundaryLayerColumnInput& input);

}  // namespace mps
