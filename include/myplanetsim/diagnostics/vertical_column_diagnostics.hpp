#pragma once

#include <cstdint>

#include "myplanetsim/config/experiment_config.hpp"
#include "myplanetsim/dynamics/vertical_column_state.hpp"
#include "myplanetsim/vertical/hybrid_pressure_coordinate.hpp"
#include "myplanetsim/vertical/vertical_transport.hpp"

namespace mps::diagnostics {

struct VerticalColumnDiagnostics {
  Real dry_mass_kg_m2;
  Real expected_dry_mass_kg_m2;
  Real dry_mass_structural_residual_kg_m2;
  Real dry_mass_budget_residual_kg_m2;
  Real potential_temperature_mass_k_kg_m2;
  Real potential_temperature_mass_budget_residual_k_kg_m2;
  Real tracer_mass_kg_m2;
  Real tracer_mass_budget_residual_kg_m2;
  Real minimum_delta_pressure_pa;
  Real maximum_delta_pressure_pa;
  Real minimum_air_mass_kg_m2;
  Real maximum_air_mass_kg_m2;
  Real minimum_theta_k;
  Real maximum_theta_k;
  Real minimum_temperature_k;
  Real maximum_temperature_k;
  Real minimum_tracer;
  Real maximum_tracer;
  Real top_mass_flux_kg_m2_s;
  Real surface_mass_flux_kg_m2_s;
  Real continuity_residual_pa_s;
  Real maximum_cfl;
  Real hydrostatic_l1_residual_m2_s2;
  Real hydrostatic_l2_residual_m2_s2;
  Real hydrostatic_linf_residual_m2_s2;
  std::uint64_t non_finite_count;
};

[[nodiscard]] VerticalColumnDiagnostics diagnose_vertical_column(
    const ExperimentConfig& config, const AtmosphericHybridCoordinate& coordinate,
    const VerticalColumnState& state, const VerticalMassFlux& mass_flux,
    Real maximum_cfl, const VerticalColumnBudget& budget);

}  // namespace mps::diagnostics
