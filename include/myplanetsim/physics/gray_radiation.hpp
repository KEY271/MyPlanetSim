#pragma once

#include <span>
#include <vector>

#include "myplanetsim/config/experiment_config.hpp"
#include "myplanetsim/core/types.hpp"

namespace mps {

// All interface arrays use top-to-bottom half-level order, j=0..K. Directional
// fluxes are non-negative magnitudes; GrayRadiationColumn::net_flux_w_m2 is upward
// positive (ADR 0017).
struct GrayRadiationColumnInput {
  std::span<const Real> pressure_half_pa;
  std::span<const Real> temperature_k;
  Real surface_temperature_k = 0.0;
  Real gravity_m_s2 = 0.0;
  Real stellar_flux_w_m2 = 0.0;
  Real cosine_solar_zenith = 0.0;
  Real surface_albedo = 0.0;
  Real surface_emissivity = 0.0;
  RadiationParameters parameters;
};

struct GrayRadiationOpticalDepth {
  std::vector<Real> shortwave;
  std::vector<Real> longwave;
};

struct GrayRadiationColumn {
  GrayRadiationOpticalDepth optical_depth;
  std::vector<Real> shortwave_down_w_m2;
  std::vector<Real> shortwave_up_w_m2;
  std::vector<Real> longwave_down_w_m2;
  std::vector<Real> longwave_up_w_m2;
  std::vector<Real> net_flux_w_m2;
  std::vector<Real> shortwave_convergence_w_m2;
  std::vector<Real> longwave_convergence_w_m2;
  std::vector<Real> radiative_convergence_w_m2;
};

struct GrayRadiativeColumnSourceInput {
  GrayRadiationColumnInput radiation;
  std::span<const Real> exner_full;
  Real heat_capacity_cp_j_kg_k = 0.0;
  Real surface_heat_capacity_j_m2_k = 0.0;
  Real air_exchange_coefficient_w_m2_k = 0.0;
  Real internal_heat_flux_w_m2 = 0.0;
};

struct GrayRadiativeColumnTendency {
  GrayRadiationColumn fluxes;
  std::vector<Real> potential_temperature_mass_k_kg_m2_s;
  Real surface_temperature_k_s = 0.0;
  Real sensible_to_atmosphere_w_m2 = 0.0;
  Real surface_storage_rate_w_m2 = 0.0;
  Real temperature_rate_bound_s_1 = 0.0;
  Real stable_time_step_s = 0.0;
};

struct GrayRadiativeColumnWorkspace {
  std::vector<Real> downward_temperature_sensitivity_w_m2_k;
  std::vector<Real> upward_temperature_sensitivity_w_m2_k;
};

// Computes physical vertical optical depths only. Angular diffusivity factors are
// applied later by the transport sweeps and are deliberately absent here.
[[nodiscard]] GrayRadiationOpticalDepth gray_radiation_optical_depth(
    std::span<const Real> pressure_half_pa, Real gravity_m_s2,
    const RadiationParameters& parameters);
void gray_radiation_optical_depth(std::span<const Real> pressure_half_pa,
                                  Real gravity_m_s2,
                                  const RadiationParameters& parameters,
                                  GrayRadiationOpticalDepth& result);

[[nodiscard]] GrayRadiationColumn gray_radiation_column(
    const GrayRadiationColumnInput& input);
void gray_radiation_column(const GrayRadiationColumnInput& input,
                           GrayRadiationColumn& result);

[[nodiscard]] GrayRadiativeColumnTendency gray_radiative_column_tendency(
    const GrayRadiativeColumnSourceInput& input);
void gray_radiative_column_tendency(const GrayRadiativeColumnSourceInput& input,
                                    GrayRadiativeColumnTendency& result);
void gray_radiative_column_tendency(const GrayRadiativeColumnSourceInput& input,
                                    GrayRadiativeColumnTendency& result,
                                    GrayRadiativeColumnWorkspace& workspace);

}  // namespace mps
