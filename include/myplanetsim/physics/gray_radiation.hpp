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
  std::vector<Real> radiative_convergence_w_m2;
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

}  // namespace mps
