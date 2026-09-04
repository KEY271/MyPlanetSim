#include "myplanetsim/diagnostics/vertical_column_diagnostics.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>

#include "myplanetsim/vertical/hydrostatic_column.hpp"

namespace mps::diagnostics {

VerticalColumnDiagnostics diagnose_vertical_column(
    const ExperimentConfig& config, const AtmosphericHybridCoordinate& coordinate,
    const VerticalColumnState& state) {
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
  return {.dry_mass_kg_m2 = std::accumulate(geometry.air_mass_kg_m2.begin(),
                                            geometry.air_mass_kg_m2.end(), 0.0),
          .expected_dry_mass_kg_m2 =
              (state.surface_pressure_pa - coordinate.top_pressure_pa()) /
              config.planet.gravity_m_s2,
          .potential_temperature_mass_k_kg_m2 =
              std::accumulate(state.potential_temperature_mass_k_kg_m2.begin(),
                              state.potential_temperature_mass_k_kg_m2.end(), 0.0),
          .tracer_mass_kg_m2 = std::accumulate(state.tracer_mass_kg_m2.begin(),
                                               state.tracer_mass_kg_m2.end(), 0.0),
          .minimum_delta_pressure_pa = *std::min_element(
              geometry.delta_pressure_pa.begin(), geometry.delta_pressure_pa.end()),
          .maximum_delta_pressure_pa = *std::max_element(
              geometry.delta_pressure_pa.begin(), geometry.delta_pressure_pa.end()),
          .minimum_air_mass_kg_m2 = *std::min_element(geometry.air_mass_kg_m2.begin(),
                                                      geometry.air_mass_kg_m2.end()),
          .minimum_theta_k = *std::min_element(theta.begin(), theta.end()),
          .minimum_temperature_k =
              *std::min_element(temperature.begin(), temperature.end()),
          .maximum_temperature_k =
              *std::max_element(temperature.begin(), temperature.end()),
          .minimum_tracer = *std::min_element(tracer.begin(), tracer.end()),
          .maximum_tracer = *std::max_element(tracer.begin(), tracer.end()),
          .hydrostatic_linf_residual_m2_s2 = [&hydrostatic] {
            Real maximum = 0.0;
            for (const Real value : hydrostatic.residual_m2_s2)
              maximum = std::max(maximum, std::abs(value));
            return maximum;
          }()};
}

}  // namespace mps::diagnostics
