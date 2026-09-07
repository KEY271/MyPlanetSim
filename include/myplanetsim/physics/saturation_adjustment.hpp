#pragma once

#include <cstddef>
#include <span>

#include "myplanetsim/physics/moist_thermodynamics.hpp"

namespace mps {

struct SaturationAdjustmentDiagnostics {
  Real condensed_water_kg_m2 = 0.0;
  Real enthalpy_change_j_m2 = 0.0;
  Real maximum_supersaturation = 0.0;
  Real maximum_vapor_mixing_ratio = 0.0;
  Real maximum_vapor_pressure_fraction = 0.0;
  std::size_t adjusted_layer_count = 0;
  std::size_t nonlinear_iterations = 0;
};

[[nodiscard]] SaturationAdjustmentDiagnostics adjust_saturation_column(
    std::span<const Real> air_mass_kg_m2, std::span<const Real> pressure_pa,
    std::span<const Real> exner_full,
    std::span<Real> potential_temperature_mass_k_kg_m2,
    std::span<Real> vapor_mass_kg_m2, const DiluteMoistThermodynamics& thermodynamics);

}  // namespace mps
