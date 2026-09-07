#include "myplanetsim/physics/boundary_layer.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

#include "myplanetsim/core/validation.hpp"

namespace mps {
namespace {

void solve_tridiagonal(std::vector<Real>& lower, std::vector<Real>& diagonal,
                       std::vector<Real>& upper, std::vector<Real>& rhs,
                       std::vector<Real>& solution) {
  const std::size_t size = diagonal.size();
  if (size == 0 || lower.size() != size || upper.size() != size || rhs.size() != size)
    throw std::invalid_argument("tridiagonal solve shape mismatch");
  for (std::size_t row = 1; row < size; ++row) {
    if (!(diagonal[row - 1] > 0.0) || !std::isfinite(diagonal[row - 1]))
      throw std::runtime_error("boundary-layer matrix pivot is invalid");
    const Real factor = lower[row] / diagonal[row - 1];
    diagonal[row] -= factor * upper[row - 1];
    rhs[row] -= factor * rhs[row - 1];
  }
  if (!(diagonal.back() > 0.0) || !std::isfinite(diagonal.back()))
    throw std::runtime_error("boundary-layer matrix pivot is invalid");
  solution.resize(size);
  solution.back() = rhs.back() / diagonal.back();
  for (std::size_t reverse = size - 1; reverse > 0; --reverse) {
    const std::size_t row = reverse - 1;
    solution[row] = (rhs[row] - upper[row] * solution[row + 1]) / diagonal[row];
  }
  if (!std::ranges::all_of(solution,
                           [](const Real value) { return std::isfinite(value); }))
    throw std::runtime_error("boundary-layer solution is not finite");
}

void build_closed_matrix(const std::span<const Real> capacity,
                         const std::span<const Real> conductance,
                         const Real time_step_s, std::vector<Real>& lower,
                         std::vector<Real>& diagonal, std::vector<Real>& upper) {
  const std::size_t size = capacity.size();
  if (conductance.size() != size + 1)
    throw std::invalid_argument("boundary-layer conductance shape mismatch");
  lower.assign(size, 0.0);
  diagonal.resize(size);
  upper.assign(size, 0.0);
  for (std::size_t row = 0; row < size; ++row) {
    const Real upper_conductance = conductance[row];
    const Real lower_conductance = conductance[row + 1];
    diagonal[row] =
        capacity[row] + time_step_s * (upper_conductance + lower_conductance);
    if (row > 0) lower[row] = -time_step_s * upper_conductance;
    if (row + 1 < size) upper[row] = -time_step_s * lower_conductance;
  }
}

[[nodiscard]] Real component(const Vec3 value, const std::size_t index) {
  if (index == 0) return value.x;
  if (index == 1) return value.y;
  return value.z;
}

void set_component(Vec3& value, const std::size_t index, const Real component_value) {
  if (index == 0)
    value.x = component_value;
  else if (index == 1)
    value.y = component_value;
  else
    value.z = component_value;
}

}  // namespace

void implicit_boundary_layer_column(const BoundaryLayerColumnInput& input,
                                    BoundaryLayerColumnResult& result,
                                    BoundaryLayerColumnWorkspace& workspace) {
  const std::size_t levels = input.potential_temperature_k.size();
  if (levels == 0 || input.velocity_m_s.size() != levels ||
      input.tracer_mixing_ratio.size() != levels ||
      input.air_mass_kg_m2.size() != levels || input.exner_full.size() != levels ||
      input.exner_half.size() != levels + 1 || input.height_full_m.size() != levels ||
      input.density_half_kg_m3.size() != levels + 1 ||
      input.eddy_diffusivity_momentum_m2_s.size() != levels + 1 ||
      input.eddy_diffusivity_heat_m2_s.size() != levels + 1 ||
      input.eddy_diffusivity_tracer_m2_s.size() != levels + 1)
    throw std::invalid_argument("boundary-layer column shape mismatch");
  require_positive(input.surface_temperature_k, "surface temperature");
  require_positive(input.surface_exner, "surface Exner function");
  require_positive(input.surface_heat_capacity_j_m2_k, "surface heat capacity");
  require_non_negative(input.surface_heat_conductance_w_m2_k,
                       "surface heat conductance");
  require_non_negative(input.surface_drag_conductance_kg_m2_s,
                       "surface drag conductance");
  require_positive(input.heat_capacity_cp_j_kg_k, "heat capacity cp");
  require_positive(input.time_step_s, "boundary-layer time step");
  if (input.eddy_diffusivity_momentum_m2_s.front() != 0.0 ||
      input.eddy_diffusivity_momentum_m2_s.back() != 0.0 ||
      input.eddy_diffusivity_heat_m2_s.front() != 0.0 ||
      input.eddy_diffusivity_heat_m2_s.back() != 0.0 ||
      input.eddy_diffusivity_tracer_m2_s.front() != 0.0 ||
      input.eddy_diffusivity_tracer_m2_s.back() != 0.0)
    throw std::invalid_argument("boundary-layer endpoint diffusivities must be zero");

  workspace.momentum_conductance.assign(levels + 1, 0.0);
  workspace.heat_conductance.assign(levels + 2, 0.0);
  workspace.tracer_conductance.assign(levels + 1, 0.0);
  workspace.momentum_capacity.resize(levels);
  workspace.heat_capacity.resize(levels + 1);
  workspace.tracer_capacity.resize(levels);
  for (std::size_t level = 0; level < levels; ++level) {
    require_positive(input.potential_temperature_k[level], "potential temperature");
    require_positive(input.air_mass_kg_m2[level], "air mass");
    require_positive(input.exner_full[level], "full-level Exner function");
    require_positive(input.exner_half[level], "half-level Exner function");
    require_positive(input.density_half_kg_m3[level], "half-level density");
    require_finite(input.height_full_m[level], "full-level height");
    require_finite(input.tracer_mixing_ratio[level], "tracer mixing ratio");
    if (!is_finite(input.velocity_m_s[level]))
      throw std::invalid_argument("velocity must be finite");
    require_non_negative(input.eddy_diffusivity_momentum_m2_s[level],
                         "momentum diffusivity");
    require_non_negative(input.eddy_diffusivity_heat_m2_s[level], "heat diffusivity");
    require_non_negative(input.eddy_diffusivity_tracer_m2_s[level],
                         "tracer diffusivity");
    workspace.momentum_capacity[level] = input.air_mass_kg_m2[level];
    workspace.tracer_capacity[level] = input.air_mass_kg_m2[level];
    workspace.heat_capacity[level] = input.heat_capacity_cp_j_kg_k *
                                     input.air_mass_kg_m2[level] *
                                     input.exner_full[level];
    if (level + 1 < levels) {
      const std::size_t interface = level + 1;
      const Real distance = input.height_full_m[level] - input.height_full_m[level + 1];
      require_positive(distance, "vertical level spacing");
      workspace.momentum_conductance[interface] =
          input.density_half_kg_m3[interface] *
          input.eddy_diffusivity_momentum_m2_s[interface] / distance;
      workspace.heat_conductance[interface] =
          input.density_half_kg_m3[interface] * input.heat_capacity_cp_j_kg_k *
          input.exner_half[interface] * input.eddy_diffusivity_heat_m2_s[interface] /
          distance;
      workspace.tracer_conductance[interface] =
          input.density_half_kg_m3[interface] *
          input.eddy_diffusivity_tracer_m2_s[interface] / distance;
    }
  }
  require_positive(input.exner_half.back(), "half-level Exner function");
  require_positive(input.density_half_kg_m3.back(), "half-level density");
  require_non_negative(input.eddy_diffusivity_momentum_m2_s.back(),
                       "momentum diffusivity");
  require_non_negative(input.eddy_diffusivity_heat_m2_s.back(), "heat diffusivity");
  require_non_negative(input.eddy_diffusivity_tracer_m2_s.back(), "tracer diffusivity");

  result.velocity_m_s.assign(input.velocity_m_s.begin(), input.velocity_m_s.end());
  workspace.closed_momentum_conductance = workspace.momentum_conductance;
  workspace.closed_momentum_conductance.back() = input.surface_drag_conductance_kg_m2_s;
  for (std::size_t axis = 0; axis < 3; ++axis) {
    build_closed_matrix(workspace.momentum_capacity,
                        workspace.closed_momentum_conductance, input.time_step_s,
                        workspace.lower, workspace.diagonal, workspace.upper);
    workspace.rhs.resize(levels);
    for (std::size_t level = 0; level < levels; ++level)
      workspace.rhs[level] = workspace.momentum_capacity[level] *
                             component(input.velocity_m_s[level], axis);
    solve_tridiagonal(workspace.lower, workspace.diagonal, workspace.upper,
                      workspace.rhs, workspace.solution);
    for (std::size_t level = 0; level < levels; ++level)
      set_component(result.velocity_m_s[level], axis, workspace.solution[level]);
  }

  result.dissipated_heat_j_m2.assign(levels, 0.0);
  Real initial_kinetic_energy = 0.0;
  Real final_kinetic_energy = 0.0;
  Real physical_shear = 0.0;
  for (std::size_t level = 0; level < levels; ++level) {
    initial_kinetic_energy +=
        0.5 * input.air_mass_kg_m2[level] * norm_squared(input.velocity_m_s[level]);
    final_kinetic_energy +=
        0.5 * input.air_mass_kg_m2[level] * norm_squared(result.velocity_m_s[level]);
    const Real numerical =
        0.5 * input.air_mass_kg_m2[level] *
        norm_squared(result.velocity_m_s[level] - input.velocity_m_s[level]);
    result.dissipated_heat_j_m2[level] += numerical;
  }
  for (std::size_t interface = 1; interface < levels; ++interface) {
    const Real dissipated = input.time_step_s *
                            workspace.momentum_conductance[interface] *
                            norm_squared(result.velocity_m_s[interface - 1] -
                                         result.velocity_m_s[interface]);
    physical_shear += dissipated;
    result.dissipated_heat_j_m2[interface - 1] += 0.5 * dissipated;
    result.dissipated_heat_j_m2[interface] += 0.5 * dissipated;
  }
  const Real physical_surface = input.time_step_s *
                                input.surface_drag_conductance_kg_m2_s *
                                norm_squared(result.velocity_m_s.back());
  result.dissipated_heat_j_m2.back() += physical_surface;
  Real numerical_dissipation = 0.0;
  for (std::size_t level = 0; level < levels; ++level) {
    numerical_dissipation +=
        0.5 * input.air_mass_kg_m2[level] *
        norm_squared(result.velocity_m_s[level] - input.velocity_m_s[level]);
  }
  const Real returned_heat = physical_shear + physical_surface + numerical_dissipation;

  result.tracer_mixing_ratio.assign(input.tracer_mixing_ratio.begin(),
                                    input.tracer_mixing_ratio.end());
  build_closed_matrix(workspace.tracer_capacity, workspace.tracer_conductance,
                      input.time_step_s, workspace.lower, workspace.diagonal,
                      workspace.upper);
  workspace.rhs.resize(levels);
  for (std::size_t level = 0; level < levels; ++level)
    workspace.rhs[level] =
        workspace.tracer_capacity[level] * input.tracer_mixing_ratio[level];
  solve_tridiagonal(workspace.lower, workspace.diagonal, workspace.upper, workspace.rhs,
                    workspace.solution);
  result.tracer_mixing_ratio = workspace.solution;
  Real initial_tracer_mass = 0.0;
  Real solved_tracer_mass = 0.0;
  Real column_mass = 0.0;
  for (std::size_t level = 0; level < levels; ++level) {
    initial_tracer_mass +=
        input.air_mass_kg_m2[level] * input.tracer_mixing_ratio[level];
    solved_tracer_mass +=
        input.air_mass_kg_m2[level] * result.tracer_mixing_ratio[level];
    column_mass += input.air_mass_kg_m2[level];
  }
  // The closed matrix conserves this mode algebraically. Project the Thomas-solve
  // roundoff back onto that null mode so very stiff columns retain the invariant.
  const Real tracer_roundoff_correction =
      (initial_tracer_mass - solved_tracer_mass) / column_mass;
  for (auto& value : result.tracer_mixing_ratio) value += tracer_roundoff_correction;

  workspace.heat_capacity.back() =
      input.surface_heat_capacity_j_m2_k * input.surface_exner;
  workspace.heat_conductance[levels] = input.surface_heat_conductance_w_m2_k;
  build_closed_matrix(workspace.heat_capacity, workspace.heat_conductance,
                      input.time_step_s, workspace.lower, workspace.diagonal,
                      workspace.upper);
  workspace.rhs.resize(levels + 1);
  for (std::size_t level = 0; level < levels; ++level)
    workspace.rhs[level] =
        workspace.heat_capacity[level] * input.potential_temperature_k[level] +
        result.dissipated_heat_j_m2[level];
  workspace.rhs.back() =
      input.surface_heat_capacity_j_m2_k * input.surface_temperature_k;
  solve_tridiagonal(workspace.lower, workspace.diagonal, workspace.upper, workspace.rhs,
                    workspace.solution);
  result.potential_temperature_k.assign(
      workspace.solution.begin(),
      workspace.solution.begin() + static_cast<std::ptrdiff_t>(levels));
  result.surface_temperature_k = input.surface_exner * workspace.solution.back();

  result.heat_flux_w_m2.assign(levels + 1, 0.0);
  result.momentum_flux_kg_m_s2.assign(levels + 1, {});
  result.tracer_flux_kg_m2_s.assign(levels + 1, 0.0);
  for (std::size_t interface = 1; interface < levels; ++interface) {
    result.heat_flux_w_m2[interface] = workspace.heat_conductance[interface] *
                                       (result.potential_temperature_k[interface] -
                                        result.potential_temperature_k[interface - 1]);
    result.momentum_flux_kg_m_s2[interface] =
        workspace.momentum_conductance[interface] *
        (result.velocity_m_s[interface] - result.velocity_m_s[interface - 1]);
    result.tracer_flux_kg_m2_s[interface] = workspace.tracer_conductance[interface] *
                                            (result.tracer_mixing_ratio[interface] -
                                             result.tracer_mixing_ratio[interface - 1]);
  }
  result.heat_flux_w_m2.back() = input.surface_heat_conductance_w_m2_k *
                                 (result.surface_temperature_k / input.surface_exner -
                                  result.potential_temperature_k.back());
  result.momentum_flux_kg_m_s2.back() =
      -input.surface_drag_conductance_kg_m2_s * result.velocity_m_s.back();

  Real atmospheric_heat_change = 0.0;
  Vec3 momentum_change{};
  Real tracer_change = 0.0;
  for (std::size_t level = 0; level < levels; ++level) {
    atmospheric_heat_change +=
        workspace.heat_capacity[level] *
        (result.potential_temperature_k[level] - input.potential_temperature_k[level]);
    momentum_change =
        momentum_change + input.air_mass_kg_m2[level] *
                              (result.velocity_m_s[level] - input.velocity_m_s[level]);
    tracer_change += input.air_mass_kg_m2[level] * (result.tracer_mixing_ratio[level] -
                                                    input.tracer_mixing_ratio[level]);
  }
  const Real surface_heat_change =
      input.surface_heat_capacity_j_m2_k *
      (result.surface_temperature_k - input.surface_temperature_k);
  const Vec3 surface_impulse = input.time_step_s * result.momentum_flux_kg_m_s2.back();
  const Real kinetic_change = final_kinetic_energy - initial_kinetic_energy;
  result.diagnostics = {
      .atmospheric_heat_change_j_m2 = atmospheric_heat_change,
      .surface_heat_change_j_m2 = surface_heat_change,
      .physical_shear_dissipation_j_m2 = physical_shear,
      .physical_surface_drag_dissipation_j_m2 = physical_surface,
      .backward_euler_dissipation_j_m2 = numerical_dissipation,
      .returned_dissipation_heat_j_m2 = returned_heat,
      .heat_budget_residual_j_m2 =
          atmospheric_heat_change + surface_heat_change - returned_heat,
      .kinetic_energy_change_j_m2 = kinetic_change,
      .kinetic_energy_identity_residual_j_m2 = -kinetic_change - returned_heat,
      .atmospheric_momentum_change_kg_m_s = momentum_change,
      .surface_stress_impulse_kg_m_s = surface_impulse,
      .momentum_budget_residual_kg_m_s = momentum_change - surface_impulse,
      .tracer_mass_change_kg_m2 = tracer_change,
  };
}

BoundaryLayerColumnResult implicit_boundary_layer_column(
    const BoundaryLayerColumnInput& input) {
  BoundaryLayerColumnResult result;
  BoundaryLayerColumnWorkspace workspace;
  implicit_boundary_layer_column(input, result, workspace);
  return result;
}

}  // namespace mps
