#pragma once

#include <cstddef>
#include <span>
#include <vector>

#include "myplanetsim/physics/moist_thermodynamics.hpp"

namespace mps {

enum class MoistConvectionBranch { kNone, kDeep, kShallow };

enum class MoistConvectionReason {
  kAdjusted,
  kNoCape,
  kThermalCondition,
  kNoShallowSolution,
};

struct SimpleBettsMillerInput {
  std::span<const Real> temperature_k;
  std::span<const Real> vapor_mixing_ratio;
  std::span<const Real> air_mass_kg_m2;
  std::span<const Real> pressure_full_pa;
  std::span<const Real> pressure_half_pa;
  Real gravity_m_s2 = 0.0;
  Real relative_humidity_reference = 0.8;
  Real relaxation_time_s = 7200.0;
  Real time_step_s = 0.0;
  Real minimum_temperature_k = 150.0;
  DiluteMoistThermodynamics thermodynamics;
};

struct SimpleBettsMillerDiagnostics {
  MoistConvectionBranch branch = MoistConvectionBranch::kNone;
  MoistConvectionReason reason = MoistConvectionReason::kNoCape;
  Real cape_j_kg = 0.0;
  Real cin_j_kg = 0.0;
  Real lcl_pressure_pa = 0.0;
  Real convection_top_pressure_pa = 0.0;
  Real raw_water_excess_kg_m2 = 0.0;
  Real raw_thermal_deficit_kg_m2 = 0.0;
  Real convective_rain_kg_m2 = 0.0;
  Real column_water_change_kg_m2 = 0.0;
  Real moist_enthalpy_change_j_m2 = 0.0;
  Real maximum_temperature_increment_k = 0.0;
  bool reaches_model_top = false;
};

struct SimpleBettsMillerResult {
  std::vector<Real> temperature_k;
  std::vector<Real> vapor_mixing_ratio;
  std::vector<Real> parcel_temperature_k;
  std::vector<Real> reference_temperature_k;
  std::vector<Real> reference_vapor_mixing_ratio;
  std::vector<Real> participation_fraction;
  SimpleBettsMillerDiagnostics diagnostics;
};

void simple_betts_miller_adjustment(const SimpleBettsMillerInput& input,
                                    SimpleBettsMillerResult& result);

[[nodiscard]] SimpleBettsMillerResult simple_betts_miller_adjustment(
    const SimpleBettsMillerInput& input);

}  // namespace mps
