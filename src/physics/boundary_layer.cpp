#include "myplanetsim/physics/boundary_layer.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "myplanetsim/core/validation.hpp"
#include "myplanetsim/physics/surface_hydrology.hpp"

namespace mps {
namespace {

constexpr Real kVonKarman = 0.4;

[[nodiscard]] Real richardson_number(const Real gravity_m_s2, const Real height_m,
                                     const Real theta_difference_k,
                                     const Real theta_reference_k,
                                     const Real speed_squared_m2_s2) {
  if (speed_squared_m2_s2 == 0.0) {
    if (theta_difference_k > 0.0) return std::numeric_limits<Real>::infinity();
    if (theta_difference_k < 0.0) return -std::numeric_limits<Real>::infinity();
    return 0.0;
  }
  return gravity_m_s2 * height_m * theta_difference_k /
         (theta_reference_k * speed_squared_m2_s2);
}

[[nodiscard]] Real stability_factor(const Real richardson,
                                    const Real critical_richardson) {
  if (richardson <= 0.0) return 1.0;
  if (richardson >= critical_richardson) return 0.0;
  const Real remaining = 1.0 - richardson / critical_richardson;
  return remaining * remaining;
}

struct TileTransfer {
  Real drag_coefficient = 0.0;
  Real heat_coefficient = 0.0;
  Real friction_velocity_m_s = 0.0;
};

[[nodiscard]] TileTransfer tile_transfer(const SurfaceRoughness roughness,
                                         const Real measurement_height_m,
                                         const Real stability,
                                         const Real effective_speed_m_s) {
  require_positive(roughness.momentum_m, "momentum roughness");
  require_positive(roughness.heat_m, "heat roughness");
  if (!(roughness.momentum_m < measurement_height_m) ||
      !(roughness.heat_m < measurement_height_m))
    throw std::invalid_argument("surface roughness must be below the lowest level");
  const Real momentum_log = std::log(measurement_height_m / roughness.momentum_m);
  const Real heat_log = std::log(measurement_height_m / roughness.heat_m);
  const Real drag = kVonKarman * kVonKarman * stability / (momentum_log * momentum_log);
  const Real heat = kVonKarman * kVonKarman * stability / (momentum_log * heat_log);
  return {.drag_coefficient = drag,
          .heat_coefficient = heat,
          .friction_velocity_m_s = std::sqrt(drag) * effective_speed_m_s};
}

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

void diagnose_boundary_layer_column(const BoundaryLayerBulkInput& input,
                                    BoundaryLayerBulkResult& result) {
  const std::size_t levels = input.potential_temperature_k.size();
  if (levels == 0 || input.temperature_k.size() != levels ||
      input.velocity_m_s.size() != levels ||
      input.pressure_half_pa.size() != levels + 1 ||
      input.exner_half.size() != levels + 1 ||
      input.height_half_m.size() != levels + 1 || input.height_full_m.size() != levels)
    throw std::invalid_argument("boundary-layer bulk shape mismatch");
  require_positive(input.surface_temperature_k, "surface temperature");
  require_positive(input.surface_exner, "surface Exner function");
  require_positive(input.gravity_m_s2, "gravity");
  require_positive(input.gas_constant_j_kg_k, "gas constant");
  require_positive(input.heat_capacity_cp_j_kg_k, "heat capacity cp");
  require_positive(input.critical_richardson, "critical Richardson number");
  require_positive(input.turbulent_prandtl, "turbulent Prandtl number");
  require_non_negative(input.gustiness_m_s, "gustiness");
  require_finite(input.land_fraction, "land fraction");
  if (input.land_fraction < 0.0 || input.land_fraction > 1.0)
    throw std::invalid_argument("land fraction must be in [0, 1]");

  for (std::size_t interface = 0; interface <= levels; ++interface) {
    require_positive(input.pressure_half_pa[interface], "half-level pressure");
    require_positive(input.exner_half[interface], "half-level Exner function");
    require_finite(input.height_half_m[interface], "half-level height");
    if (interface > 0) {
      if (!(input.pressure_half_pa[interface] > input.pressure_half_pa[interface - 1]))
        throw std::invalid_argument("half-level pressure must increase downward");
      if (!(input.height_half_m[interface] < input.height_half_m[interface - 1]))
        throw std::invalid_argument("half-level height must decrease downward");
    }
  }
  if (input.height_half_m.back() != 0.0)
    throw std::invalid_argument("boundary-layer surface height must be zero");
  for (std::size_t level = 0; level < levels; ++level) {
    require_positive(input.potential_temperature_k[level], "potential temperature");
    require_positive(input.temperature_k[level], "temperature");
    require_positive(input.height_full_m[level], "full-level height");
    if (!is_finite(input.velocity_m_s[level]))
      throw std::invalid_argument("velocity must be finite");
    if (level > 0 && !(input.height_full_m[level] < input.height_full_m[level - 1]))
      throw std::invalid_argument("full-level height must decrease downward");
  }

  const std::size_t bottom = levels - 1;
  const Real theta_surface = input.surface_temperature_k / input.surface_exner;
  const Real theta_air = input.potential_temperature_k[bottom];
  const Real theta_reference = 0.5 * (theta_surface + theta_air);
  require_positive(theta_reference, "surface reference potential temperature");
  const Real wind_squared = norm_squared(input.velocity_m_s[bottom]);
  const Real effective_speed_squared =
      wind_squared + input.gustiness_m_s * input.gustiness_m_s;
  const Real effective_speed = std::sqrt(effective_speed_squared);
  const Real surface_richardson = richardson_number(
      input.gravity_m_s2, input.height_full_m[bottom], theta_air - theta_surface,
      theta_reference, effective_speed_squared);
  const Real surface_stability =
      stability_factor(surface_richardson, input.critical_richardson);
  const auto land = tile_transfer(input.land_roughness, input.height_full_m[bottom],
                                  surface_stability, effective_speed);
  const auto ocean = tile_transfer(input.ocean_roughness, input.height_full_m[bottom],
                                   surface_stability, effective_speed);
  const Real ocean_fraction = 1.0 - input.land_fraction;
  const Real drag_coefficient = input.land_fraction * land.drag_coefficient +
                                ocean_fraction * ocean.drag_coefficient;
  const Real heat_coefficient = input.land_fraction * land.heat_coefficient +
                                ocean_fraction * ocean.heat_coefficient;

  bool shallow_unresolved = false;
  bool reaches_top = false;
  Real boundary_layer_height = input.height_full_m[bottom];
  if (surface_richardson >= input.critical_richardson) {
    shallow_unresolved = true;
  } else {
    Real previous_height = input.height_full_m[bottom];
    Real previous_richardson = 0.0;
    bool crossed = false;
    for (std::size_t reverse = bottom; reverse > 0; --reverse) {
      const std::size_t level = reverse - 1;
      const Real reference = 0.5 * (input.potential_temperature_k[level] + theta_air);
      const Real speed_squared = norm_squared(input.velocity_m_s[level]) +
                                 input.gustiness_m_s * input.gustiness_m_s;
      const Real bulk_richardson = richardson_number(
          input.gravity_m_s2, input.height_full_m[level],
          input.potential_temperature_k[level] - theta_air, reference, speed_squared);
      if (bulk_richardson >= input.critical_richardson) {
        Real fraction = 0.0;
        if (std::isfinite(bulk_richardson))
          fraction = std::clamp((input.critical_richardson - previous_richardson) /
                                    (bulk_richardson - previous_richardson),
                                0.0, 1.0);
        boundary_layer_height =
            previous_height + fraction * (input.height_full_m[level] - previous_height);
        crossed = true;
        break;
      }
      previous_height = input.height_full_m[level];
      previous_richardson = bulk_richardson;
    }
    if (!crossed) {
      boundary_layer_height = input.height_half_m.front();
      reaches_top = true;
    }
  }

  result.density_half_kg_m3.resize(levels + 1);
  result.eddy_diffusivity_momentum_m2_s.assign(levels + 1, 0.0);
  result.eddy_diffusivity_heat_m2_s.assign(levels + 1, 0.0);
  result.eddy_diffusivity_tracer_m2_s.assign(levels + 1, 0.0);
  for (std::size_t interface = 0; interface <= levels; ++interface) {
    Real interface_temperature = input.surface_temperature_k;
    if (interface == 0)
      interface_temperature = input.temperature_k.front();
    else if (interface < levels)
      interface_temperature =
          0.5 * (input.temperature_k[interface - 1] + input.temperature_k[interface]);
    result.density_half_kg_m3[interface] =
        input.pressure_half_pa[interface] /
        (input.gas_constant_j_kg_k * interface_temperature);
    if (interface == 0 || interface == levels) continue;
    const Real height = input.height_half_m[interface];
    if (!(height > 0.0) || !(height < boundary_layer_height)) continue;
    const Real shape = height * (1.0 - height / boundary_layer_height) *
                       (1.0 - height / boundary_layer_height);
    const Real land_k = kVonKarman * land.friction_velocity_m_s * shape;
    const Real ocean_k = kVonKarman * ocean.friction_velocity_m_s * shape;
    const Real momentum_k = input.land_fraction * land_k + ocean_fraction * ocean_k;
    result.eddy_diffusivity_momentum_m2_s[interface] = momentum_k;
    result.eddy_diffusivity_heat_m2_s[interface] = momentum_k / input.turbulent_prandtl;
    result.eddy_diffusivity_tracer_m2_s[interface] =
        result.eddy_diffusivity_heat_m2_s[interface];
  }

  const Real surface_density = result.density_half_kg_m3.back();
  result.surface_heat_conductance_w_m2_k =
      surface_density * input.heat_capacity_cp_j_kg_k * heat_coefficient *
      effective_speed * input.surface_exner;
  result.surface_drag_conductance_kg_m2_s =
      surface_density * drag_coefficient * effective_speed;
  result.surface_water_conductance_land_kg_m2_s =
      surface_density * land.heat_coefficient * effective_speed;
  result.surface_water_conductance_ocean_kg_m2_s =
      surface_density * ocean.heat_coefficient * effective_speed;
  const Real sensible_heat =
      result.surface_heat_conductance_w_m2_k * (theta_surface - theta_air);
  const Vec3 stress =
      -result.surface_drag_conductance_kg_m2_s * input.velocity_m_s[bottom];
  result.diagnostics = {
      .surface_richardson = surface_richardson,
      .surface_stability_factor = surface_stability,
      .drag_coefficient = drag_coefficient,
      .heat_exchange_coefficient = heat_coefficient,
      .boundary_layer_height_m = boundary_layer_height,
      .maximum_momentum_diffusivity_m2_s =
          *std::max_element(result.eddy_diffusivity_momentum_m2_s.begin(),
                            result.eddy_diffusivity_momentum_m2_s.end()),
      .maximum_heat_diffusivity_m2_s =
          *std::max_element(result.eddy_diffusivity_heat_m2_s.begin(),
                            result.eddy_diffusivity_heat_m2_s.end()),
      .sensible_heat_flux_w_m2 = sensible_heat,
      .surface_stress_kg_m_s2 = stress,
      .shallow_stable_layer_unresolved = shallow_unresolved,
      .reaches_model_top = reaches_top,
  };
}

BoundaryLayerBulkResult diagnose_boundary_layer_column(
    const BoundaryLayerBulkInput& input) {
  BoundaryLayerBulkResult result;
  diagnose_boundary_layer_column(input, result);
  return result;
}

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
  if (input.enable_surface_water_exchange) {
    input.moist_thermodynamics.validate();
    require_positive(input.surface_pressure_pa, "surface pressure");
    require_finite(input.land_fraction, "land fraction");
    if (input.land_fraction < 0.0 || input.land_fraction > 1.0)
      throw std::invalid_argument("land fraction must be in [0, 1]");
    require_non_negative(input.land_water_kg_m2, "land water");
    require_non_negative(input.bucket_capacity_kg_m2, "bucket capacity");
    require_non_negative(input.surface_water_conductance_land_kg_m2_s,
                         "land water conductance");
    require_non_negative(input.surface_water_conductance_ocean_kg_m2_s,
                         "ocean water conductance");
  }
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
  const Real old_surface_theta = input.surface_temperature_k / input.surface_exner;
  for (std::size_t level = 0; level < levels; ++level) {
    Real old_flux_convergence = 0.0;
    if (level > 0)
      old_flux_convergence += workspace.heat_conductance[level] *
                              (input.potential_temperature_k[level - 1] -
                               input.potential_temperature_k[level]);
    if (level + 1 < levels)
      old_flux_convergence += workspace.heat_conductance[level + 1] *
                              (input.potential_temperature_k[level + 1] -
                               input.potential_temperature_k[level]);
    else
      old_flux_convergence +=
          workspace.heat_conductance[levels] *
          (old_surface_theta - input.potential_temperature_k[level]);
    workspace.rhs[level] =
        input.time_step_s * old_flux_convergence + result.dissipated_heat_j_m2[level];
  }
  workspace.rhs.back() = input.time_step_s * workspace.heat_conductance[levels] *
                         (input.potential_temperature_k.back() - old_surface_theta);
  solve_tridiagonal(workspace.lower, workspace.diagonal, workspace.upper, workspace.rhs,
                    workspace.solution);
  workspace.heat_increment = workspace.solution;
  result.potential_temperature_k.resize(levels);
  for (std::size_t level = 0; level < levels; ++level)
    result.potential_temperature_k[level] =
        input.potential_temperature_k[level] + workspace.solution[level];
  result.surface_temperature_k =
      input.surface_temperature_k + input.surface_exner * workspace.solution.back();

  SurfaceEvaporationPartition surface_water;
  std::size_t surface_water_iterations = 0;
  if (input.enable_surface_water_exchange) {
    workspace.tracer_surface_response.resize(levels);
    build_closed_matrix(workspace.tracer_capacity, workspace.tracer_conductance,
                        input.time_step_s, workspace.lower, workspace.diagonal,
                        workspace.upper);
    workspace.rhs.assign(levels, 0.0);
    workspace.rhs.back() = input.time_step_s;
    solve_tridiagonal(workspace.lower, workspace.diagonal, workspace.upper,
                      workspace.rhs, workspace.tracer_surface_response);

    workspace.heat_surface_response.resize(levels + 1);
    build_closed_matrix(workspace.heat_capacity, workspace.heat_conductance,
                        input.time_step_s, workspace.lower, workspace.diagonal,
                        workspace.upper);
    workspace.rhs.assign(levels + 1, 0.0);
    workspace.rhs.back() =
        -input.time_step_s * input.moist_thermodynamics.latent_heat_vaporization_j_kg;
    solve_tridiagonal(workspace.lower, workspace.diagonal, workspace.upper,
                      workspace.rhs, workspace.heat_surface_response);

    const auto partition = [&](const Real cell_flux) {
      const Real surface_temperature =
          result.surface_temperature_k +
          input.surface_exner * workspace.heat_surface_response.back() * cell_flux;
      const Real bottom_vapor = result.tracer_mixing_ratio.back() +
                                workspace.tracer_surface_response.back() * cell_flux;
      const Real deficit =
          saturation_mixing_ratio(surface_temperature, input.surface_pressure_pa,
                                  input.moist_thermodynamics) -
          bottom_vapor;
      return partition_surface_evaporation(
          deficit, input.surface_water_conductance_land_kg_m2_s,
          input.surface_water_conductance_ocean_kg_m2_s, input.land_fraction,
          input.land_water_kg_m2, input.bucket_capacity_kg_m2,
          input.bucket_wet_threshold_fraction, input.time_step_s);
    };
    const auto residual = [&](const Real cell_flux) {
      return cell_flux - partition(cell_flux).cell_flux_kg_m2_s;
    };

    Real lower = -std::numeric_limits<Real>::infinity();
    Real upper = std::numeric_limits<Real>::infinity();
    for (std::size_t level = 0; level < levels; ++level) {
      const Real tracer_response = workspace.tracer_surface_response[level];
      if (tracer_response > 0.0)
        lower = std::max(lower,
                         -0.999 * result.tracer_mixing_ratio[level] / tracer_response);
      if (tracer_response > 0.0)
        upper =
            std::min(upper, 0.999 *
                                (input.moist_thermodynamics.maximum_vapor_mixing_ratio -
                                 result.tracer_mixing_ratio[level]) /
                                tracer_response);
      const Real heat_response = workspace.heat_surface_response[level];
      if (heat_response < 0.0)
        upper = std::min(
            upper, -0.999 * result.potential_temperature_k[level] / heat_response);
    }
    const Real surface_temperature_response =
        input.surface_exner * workspace.heat_surface_response.back();
    if (surface_temperature_response < 0.0)
      upper = std::min(
          upper, -0.999 * result.surface_temperature_k / surface_temperature_response);
    if (!std::isfinite(lower)) lower = -1.0;
    if (!std::isfinite(upper)) upper = 1.0;
    lower = std::min(lower, 0.0);
    upper = std::max(upper, 0.0);

    Real solved_flux = 0.0;
    const Real at_zero = residual(0.0);
    if (std::abs(at_zero) > 1e-15) {
      Real bracket_lower = at_zero > 0.0 ? lower : 0.0;
      Real bracket_upper = at_zero > 0.0 ? 0.0 : upper;
      Real value_lower = residual(bracket_lower);
      Real value_upper = residual(bracket_upper);
      if (!(value_lower <= 0.0 && value_upper >= 0.0))
        throw std::runtime_error("implicit surface-water flux is not bracketed");
      for (std::size_t iteration = 0; iteration < 100; ++iteration) {
        solved_flux = 0.5 * (bracket_lower + bracket_upper);
        const Real value = residual(solved_flux);
        ++surface_water_iterations;
        if (std::abs(value) <= 1e-14 + 1e-11 * std::abs(solved_flux)) break;
        if (value < 0.0)
          bracket_lower = solved_flux;
        else
          bracket_upper = solved_flux;
        if (iteration == 99)
          throw std::runtime_error("implicit surface-water flux did not converge");
      }
    }
    surface_water = partition(solved_flux);
    for (std::size_t level = 0; level < levels; ++level) {
      result.tracer_mixing_ratio[level] +=
          workspace.tracer_surface_response[level] * solved_flux;
      const Real heat_increment = workspace.heat_surface_response[level] * solved_flux;
      result.potential_temperature_k[level] += heat_increment;
      workspace.heat_increment[level] += heat_increment;
    }
    const Real surface_theta_increment =
        workspace.heat_surface_response.back() * solved_flux;
    result.surface_temperature_k += input.surface_exner * surface_theta_increment;
    workspace.heat_increment.back() += surface_theta_increment;
    result.land_water_kg_m2 = surface_water.final_land_water_kg_m2;
    result.surface_water_flux_kg_m2_s = solved_flux;
    result.land_water_flux_kg_m2_s = surface_water.land_flux_kg_m2_s;
    result.ocean_water_flux_kg_m2_s = surface_water.ocean_flux_kg_m2_s;
  } else {
    result.land_water_kg_m2 = input.land_water_kg_m2;
  }

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
  if (input.enable_surface_water_exchange)
    result.tracer_flux_kg_m2_s.back() = result.surface_water_flux_kg_m2_s;

  Real atmospheric_heat_change = 0.0;
  Vec3 momentum_change{};
  Real tracer_change = 0.0;
  for (std::size_t level = 0; level < levels; ++level) {
    atmospheric_heat_change +=
        workspace.heat_capacity[level] * workspace.heat_increment[level];
    momentum_change =
        momentum_change + input.air_mass_kg_m2[level] *
                              (result.velocity_m_s[level] - input.velocity_m_s[level]);
    tracer_change += input.air_mass_kg_m2[level] * (result.tracer_mixing_ratio[level] -
                                                    input.tracer_mixing_ratio[level]);
  }
  const Real surface_heat_change =
      workspace.heat_capacity.back() * workspace.heat_increment.back();
  const Vec3 surface_impulse = input.time_step_s * result.momentum_flux_kg_m_s2.back();
  const Real kinetic_change = final_kinetic_energy - initial_kinetic_energy;
  const Real evaporation = input.time_step_s * result.surface_water_flux_kg_m2_s;
  const Real land_evaporation = input.time_step_s * result.land_water_flux_kg_m2_s;
  const Real ocean_evaporation = input.time_step_s * result.ocean_water_flux_kg_m2_s;
  const Real cell_runoff = input.land_fraction * surface_water.runoff_kg_m2_land;
  const Real ocean_water_change =
      input.enable_surface_water_exchange && input.land_fraction < 1.0
          ? -(1.0 - input.land_fraction) * ocean_evaporation + cell_runoff
          : 0.0;
  const Real external_outflow =
      input.enable_surface_water_exchange && input.land_fraction == 1.0 ? cell_runoff
                                                                        : 0.0;
  const Real land_water_change =
      input.land_fraction * (result.land_water_kg_m2 - input.land_water_kg_m2);
  const Real latent_surface_heat =
      -input.moist_thermodynamics.latent_heat_vaporization_j_kg * evaporation;
  result.diagnostics = {
      .atmospheric_heat_change_j_m2 = atmospheric_heat_change,
      .surface_heat_change_j_m2 = surface_heat_change,
      .physical_shear_dissipation_j_m2 = physical_shear,
      .physical_surface_drag_dissipation_j_m2 = physical_surface,
      .backward_euler_dissipation_j_m2 = numerical_dissipation,
      .returned_dissipation_heat_j_m2 = returned_heat,
      .heat_budget_residual_j_m2 = atmospheric_heat_change + surface_heat_change -
                                   returned_heat - latent_surface_heat,
      .kinetic_energy_change_j_m2 = kinetic_change,
      .kinetic_energy_identity_residual_j_m2 = -kinetic_change - returned_heat,
      .atmospheric_momentum_change_kg_m_s = momentum_change,
      .surface_stress_impulse_kg_m_s = surface_impulse,
      .momentum_budget_residual_kg_m_s = momentum_change - surface_impulse,
      .tracer_mass_change_kg_m2 = tracer_change,
      .evaporation_kg_m2 = evaporation,
      .land_evaporation_kg_m2_land = land_evaporation,
      .ocean_evaporation_kg_m2_ocean = ocean_evaporation,
      .runoff_kg_m2_land = surface_water.runoff_kg_m2_land,
      .ocean_water_change_kg_m2 = ocean_water_change,
      .external_outflow_kg_m2 = external_outflow,
      .latent_surface_heat_change_j_m2 = latent_surface_heat,
      .moist_enthalpy_budget_residual_j_m2 =
          atmospheric_heat_change +
          input.moist_thermodynamics.latent_heat_vaporization_j_kg * tracer_change +
          surface_heat_change - returned_heat,
      .water_budget_residual_kg_m2 =
          tracer_change + land_water_change + ocean_water_change + external_outflow,
      .surface_water_iterations = surface_water_iterations,
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
