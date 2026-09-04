#include "myplanetsim/diagnostics/vertical_column_diagnostics.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <stdexcept>

#include "myplanetsim/vertical/hydrostatic_column.hpp"

namespace mps::diagnostics {

VerticalColumnDiagnostics diagnose_vertical_column(
    const ExperimentConfig& config, const AtmosphericHybridCoordinate& coordinate,
    const VerticalColumnState& state, const VerticalMassFlux& mass_flux,
    const Real maximum_cfl, const VerticalColumnBudget& budget) {
  if (mass_flux.interface_flux_kg_m2_s.size() != coordinate.levels() + 1 ||
      mass_flux.target_air_mass_tendency_kg_m2_s.size() != coordinate.levels()) {
    throw std::invalid_argument("vertical diagnostic mass-flux shape is invalid");
  }
  const auto geometry = coordinate.geometry(
      state.surface_pressure_pa, config.planet.gravity_m_s2,
      config.planet.gas_constant_j_kg_k, config.planet.heat_capacity_cp_j_kg_k,
      config.planet.reference_pressure_pa);
  std::vector<Real> theta(coordinate.levels());
  std::vector<Real> tracer(coordinate.levels());
  for (std::size_t k = 0; k < coordinate.levels(); ++k) {
    theta[k] = state.potential_temperature_mass_k_kg_m2[k] / geometry.air_mass_kg_m2[k];
    tracer[k] = state.tracer_mass_kg_m2[k] / geometry.air_mass_kg_m2[k];
  }
  const auto hydrostatic = integrate_hydrostatic_column(
      geometry, theta, config.planet.heat_capacity_cp_j_kg_k,
      config.planet.gravity_m_s2, config.vertical.surface_geopotential_m2_s2);
  std::vector<Real> temperature(theta.size());
  for (std::size_t k = 0; k < theta.size(); ++k) {
    temperature[k] = theta[k] * geometry.exner_full[k];
  }
  const Real dry_mass = std::accumulate(geometry.air_mass_kg_m2.begin(),
                                        geometry.air_mass_kg_m2.end(), 0.0);
  const Real expected_dry_mass =
      (state.surface_pressure_pa - coordinate.top_pressure_pa()) /
      config.planet.gravity_m_s2;
  const Real theta_mass =
      std::accumulate(state.potential_temperature_mass_k_kg_m2.begin(),
                      state.potential_temperature_mass_k_kg_m2.end(), 0.0);
  const Real tracer_mass = std::accumulate(state.tracer_mass_kg_m2.begin(),
                                           state.tracer_mass_kg_m2.end(), 0.0);
  Real hydrostatic_l1 = 0.0;
  Real hydrostatic_l2 = 0.0;
  Real hydrostatic_linf = 0.0;
  Real total_mass = 0.0;
  for (std::size_t k = 0; k < coordinate.levels(); ++k) {
    const Real residual = hydrostatic.residual_m2_s2[k];
    const Real weight = geometry.air_mass_kg_m2[k];
    hydrostatic_l1 += weight * std::abs(residual);
    hydrostatic_l2 += weight * residual * residual;
    hydrostatic_linf = std::max(hydrostatic_linf, std::abs(residual));
    total_mass += weight;
  }
  hydrostatic_l1 /= total_mass;
  hydrostatic_l2 = std::sqrt(hydrostatic_l2 / total_mass);

  std::uint64_t non_finite_count = 0;
  const auto count_non_finite = [&non_finite_count](const auto& values) {
    for (const Real value : values) {
      if (!std::isfinite(value)) {
        ++non_finite_count;
      }
    }
  };
  count_non_finite(geometry.pressure_half_pa);
  count_non_finite(geometry.exner_half);
  count_non_finite(geometry.pressure_full_pa);
  count_non_finite(geometry.exner_full);
  count_non_finite(geometry.delta_pressure_pa);
  count_non_finite(geometry.air_mass_kg_m2);
  count_non_finite(state.potential_temperature_mass_k_kg_m2);
  count_non_finite(state.tracer_mass_kg_m2);
  count_non_finite(theta);
  count_non_finite(temperature);
  count_non_finite(tracer);
  count_non_finite(hydrostatic.geopotential_half_m2_s2);
  count_non_finite(hydrostatic.geopotential_full_m2_s2);
  count_non_finite(hydrostatic.residual_m2_s2);
  count_non_finite(mass_flux.target_air_mass_tendency_kg_m2_s);
  count_non_finite(mass_flux.interface_flux_kg_m2_s);
  const auto count_non_finite_scalar = [&non_finite_count](const Real value) {
    if (!std::isfinite(value)) {
      ++non_finite_count;
    }
  };
  count_non_finite_scalar(state.time_s);
  count_non_finite_scalar(state.surface_pressure_pa);
  count_non_finite_scalar(mass_flux.surface_pressure_tendency_pa_s);
  count_non_finite_scalar(mass_flux.continuity_residual_pa_s);
  count_non_finite_scalar(maximum_cfl);
  count_non_finite_scalar(budget.initial_dry_mass_kg_m2);
  count_non_finite_scalar(budget.initial_potential_temperature_mass_k_kg_m2);
  count_non_finite_scalar(budget.initial_tracer_mass_kg_m2);
  count_non_finite_scalar(budget.integrated_dry_mass_source_kg_m2);
  count_non_finite_scalar(budget.integrated_potential_temperature_mass_source_k_kg_m2);
  count_non_finite_scalar(budget.integrated_tracer_mass_source_kg_m2);

  const Real top_mass_flux = mass_flux.interface_flux_kg_m2_s.front();
  const Real surface_mass_flux = mass_flux.interface_flux_kg_m2_s.back();
  return {.dry_mass_kg_m2 = dry_mass,
          .expected_dry_mass_kg_m2 = expected_dry_mass,
          .dry_mass_structural_residual_kg_m2 = dry_mass - expected_dry_mass,
          .dry_mass_budget_residual_kg_m2 = dry_mass - budget.initial_dry_mass_kg_m2 -
                                            budget.integrated_dry_mass_source_kg_m2,
          .potential_temperature_mass_k_kg_m2 = theta_mass,
          .potential_temperature_mass_budget_residual_k_kg_m2 =
              theta_mass - budget.initial_potential_temperature_mass_k_kg_m2 -
              budget.integrated_potential_temperature_mass_source_k_kg_m2,
          .tracer_mass_kg_m2 = tracer_mass,
          .tracer_mass_budget_residual_kg_m2 =
              tracer_mass - budget.initial_tracer_mass_kg_m2 -
              budget.integrated_tracer_mass_source_kg_m2,
          .minimum_delta_pressure_pa = *std::min_element(
              geometry.delta_pressure_pa.begin(), geometry.delta_pressure_pa.end()),
          .maximum_delta_pressure_pa = *std::max_element(
              geometry.delta_pressure_pa.begin(), geometry.delta_pressure_pa.end()),
          .minimum_air_mass_kg_m2 = *std::min_element(geometry.air_mass_kg_m2.begin(),
                                                      geometry.air_mass_kg_m2.end()),
          .maximum_air_mass_kg_m2 = *std::max_element(geometry.air_mass_kg_m2.begin(),
                                                      geometry.air_mass_kg_m2.end()),
          .minimum_theta_k = *std::min_element(theta.begin(), theta.end()),
          .maximum_theta_k = *std::max_element(theta.begin(), theta.end()),
          .minimum_temperature_k =
              *std::min_element(temperature.begin(), temperature.end()),
          .maximum_temperature_k =
              *std::max_element(temperature.begin(), temperature.end()),
          .minimum_tracer = *std::min_element(tracer.begin(), tracer.end()),
          .maximum_tracer = *std::max_element(tracer.begin(), tracer.end()),
          .top_mass_flux_kg_m2_s = top_mass_flux,
          .surface_mass_flux_kg_m2_s = surface_mass_flux,
          .continuity_residual_pa_s = mass_flux.continuity_residual_pa_s,
          .maximum_cfl = maximum_cfl,
          .hydrostatic_l1_residual_m2_s2 = hydrostatic_l1,
          .hydrostatic_l2_residual_m2_s2 = hydrostatic_l2,
          .hydrostatic_linf_residual_m2_s2 = hydrostatic_linf,
          .non_finite_count = non_finite_count};
}

}  // namespace mps::diagnostics
