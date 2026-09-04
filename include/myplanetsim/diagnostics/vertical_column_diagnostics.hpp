#pragma once

#include "myplanetsim/config/experiment_config.hpp"
#include "myplanetsim/dynamics/vertical_column_state.hpp"
#include "myplanetsim/vertical/hybrid_pressure_coordinate.hpp"

namespace mps::diagnostics {

struct VerticalColumnDiagnostics {
  Real dry_mass_kg_m2;
  Real expected_dry_mass_kg_m2;
  Real potential_temperature_mass_k_kg_m2;
  Real tracer_mass_kg_m2;
  Real minimum_delta_pressure_pa;
  Real maximum_delta_pressure_pa;
  Real minimum_air_mass_kg_m2;
  Real minimum_theta_k;
  Real minimum_temperature_k;
  Real maximum_temperature_k;
  Real minimum_tracer;
  Real maximum_tracer;
  Real hydrostatic_linf_residual_m2_s2;
};

[[nodiscard]] VerticalColumnDiagnostics diagnose_vertical_column(
    const ExperimentConfig& config, const AtmosphericHybridCoordinate& coordinate,
    const VerticalColumnState& state);

}  // namespace mps::diagnostics
