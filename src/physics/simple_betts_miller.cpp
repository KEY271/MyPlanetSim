#include "myplanetsim/physics/simple_betts_miller.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

#include "myplanetsim/core/validation.hpp"

namespace mps {
namespace {

[[nodiscard]] Real moist_adiabatic_slope(
    const Real temperature, const Real pressure,
    const DiluteMoistThermodynamics& thermodynamics) {
  const Real saturation =
      saturation_mixing_ratio(temperature, pressure, thermodynamics);
  const Real latent = thermodynamics.latent_heat_vaporization_j_kg;
  return (thermodynamics.gas_constant_dry_air_j_kg_k * temperature +
          latent * saturation) /
         (thermodynamics.heat_capacity_cp_j_kg_k +
          latent * latent * saturation /
              (thermodynamics.gas_constant_water_vapor_j_kg_k * temperature *
               temperature));
}

[[nodiscard]] Real integrate_moist_adiabat(
    const Real initial_temperature, const Real initial_pressure,
    const Real final_pressure, const DiluteMoistThermodynamics& thermodynamics) {
  const Real interval = std::log(final_pressure / initial_pressure);
  const std::size_t steps = std::max<std::size_t>(
      1, static_cast<std::size_t>(std::ceil(std::abs(interval) / 0.025)));
  const Real step = interval / static_cast<Real>(steps);
  Real temperature = initial_temperature;
  Real log_pressure = std::log(initial_pressure);
  for (std::size_t index = 0; index < steps; ++index) {
    const Real pressure = std::exp(log_pressure);
    const Real first = moist_adiabatic_slope(temperature, pressure, thermodynamics);
    const Real midpoint_temperature = temperature + 0.5 * step * first;
    const Real midpoint_pressure = std::exp(log_pressure + 0.5 * step);
    temperature += step * moist_adiabatic_slope(midpoint_temperature, midpoint_pressure,
                                                thermodynamics);
    log_pressure += step;
  }
  return temperature;
}

struct ParcelProfile {
  std::vector<Real> temperature;
  Real lcl_pressure_pa = 0.0;
};

[[nodiscard]] ParcelProfile lift_surface_parcel(const SimpleBettsMillerInput& input) {
  const std::size_t levels = input.temperature_k.size();
  const std::size_t bottom = levels - 1;
  ParcelProfile result{.temperature = std::vector<Real>(levels),
                       .lcl_pressure_pa = input.pressure_full_pa.front()};
  Real parcel_temperature = input.temperature_k[bottom];
  Real parcel_vapor = input.vapor_mixing_ratio[bottom];
  const Real bottom_saturation = saturation_mixing_ratio(
      parcel_temperature, input.pressure_full_pa[bottom], input.thermodynamics);
  if (parcel_vapor > bottom_saturation) {
    const Real warming = input.thermodynamics.latent_heat_vaporization_j_kg /
                         input.thermodynamics.heat_capacity_cp_j_kg_k;
    const auto residual = [&](const Real condensed) {
      return parcel_vapor - condensed -
             saturation_mixing_ratio(parcel_temperature + warming * condensed,
                                     input.pressure_full_pa[bottom],
                                     input.thermodynamics);
    };
    Real lower = 0.0;
    Real upper = parcel_vapor;
    for (std::size_t iteration = 0; iteration < 100; ++iteration) {
      const Real midpoint = 0.5 * (lower + upper);
      if (residual(midpoint) > 0.0)
        lower = midpoint;
      else
        upper = midpoint;
    }
    const Real condensed = 0.5 * (lower + upper);
    parcel_temperature += warming * condensed;
    parcel_vapor -= condensed;
  }
  result.temperature[bottom] = parcel_temperature;
  bool saturated =
      parcel_vapor >= saturation_mixing_ratio(parcel_temperature,
                                              input.pressure_full_pa[bottom],
                                              input.thermodynamics) -
                          1e-14;
  if (saturated) result.lcl_pressure_pa = input.pressure_full_pa[bottom];
  Real previous_pressure = input.pressure_full_pa[bottom];
  for (std::size_t reverse = bottom; reverse > 0; --reverse) {
    const std::size_t level = reverse - 1;
    const Real target_pressure = input.pressure_full_pa[level];
    if (!saturated) {
      const Real kappa = input.thermodynamics.gas_constant_dry_air_j_kg_k /
                         input.thermodynamics.heat_capacity_cp_j_kg_k;
      const Real dry_target =
          parcel_temperature * std::pow(target_pressure / previous_pressure, kappa);
      if (parcel_vapor <
          saturation_mixing_ratio(dry_target, target_pressure, input.thermodynamics)) {
        parcel_temperature = dry_target;
      } else {
        const Real base_temperature = parcel_temperature;
        const Real base_log_pressure = std::log(previous_pressure);
        Real lower_log = std::log(target_pressure);
        Real upper_log = base_log_pressure;
        for (std::size_t iteration = 0; iteration < 80; ++iteration) {
          const Real midpoint_log = 0.5 * (lower_log + upper_log);
          const Real dry_temperature =
              base_temperature * std::exp(kappa * (midpoint_log - base_log_pressure));
          const Real value = parcel_vapor - saturation_mixing_ratio(
                                                dry_temperature, std::exp(midpoint_log),
                                                input.thermodynamics);
          if (value >= 0.0)
            lower_log = midpoint_log;
          else
            upper_log = midpoint_log;
        }
        const Real lcl_log = 0.5 * (lower_log + upper_log);
        result.lcl_pressure_pa = std::exp(lcl_log);
        const Real lcl_temperature =
            base_temperature * std::exp(kappa * (lcl_log - base_log_pressure));
        parcel_temperature =
            integrate_moist_adiabat(lcl_temperature, result.lcl_pressure_pa,
                                    target_pressure, input.thermodynamics);
        saturated = true;
      }
    } else {
      parcel_temperature = integrate_moist_adiabat(
          parcel_temperature, previous_pressure, target_pressure, input.thermodynamics);
    }
    result.temperature[level] = parcel_temperature;
    previous_pressure = target_pressure;
  }
  return result;
}

void diagnose_buoyancy(const SimpleBettsMillerInput& input,
                       const std::span<const Real> parcel_temperature,
                       std::vector<Real>& participation,
                       SimpleBettsMillerDiagnostics& diagnostics) {
  const std::size_t levels = input.temperature_k.size();
  const std::size_t bottom = levels - 1;
  std::vector<Real> buoyancy(levels);
  for (std::size_t level = 0; level < levels; ++level)
    buoyancy[level] = input.gravity_m_s2 *
                      (parcel_temperature[level] - input.temperature_k[level]) /
                      input.temperature_k[level];

  bool positive_buoyancy = false;
  diagnostics.cape_j_kg = 0.0;
  diagnostics.cin_j_kg = 0.0;
  for (std::size_t reverse = bottom; reverse > 0; --reverse) {
    const std::size_t upper = reverse - 1;
    const std::size_t lower = reverse;
    const Real thickness =
        input.thermodynamics.gas_constant_dry_air_j_kg_k * 0.5 *
        (input.temperature_k[upper] + input.temperature_k[lower]) / input.gravity_m_s2 *
        std::log(input.pressure_full_pa[lower] / input.pressure_full_pa[upper]);
    const Real mean_buoyancy = 0.5 * (buoyancy[upper] + buoyancy[lower]);
    if (mean_buoyancy > 0.0) {
      diagnostics.cape_j_kg += mean_buoyancy * thickness;
      positive_buoyancy = true;
    } else {
      diagnostics.cin_j_kg -= mean_buoyancy * thickness;
    }
  }

  participation.assign(levels, 0.0);
  if (!positive_buoyancy || !(diagnostics.cape_j_kg > 1e-12)) return;
  participation.assign(levels, 1.0);
  diagnostics.reaches_model_top = true;
  diagnostics.convection_top_pressure_pa = input.pressure_half_pa.front();
  bool found_lfc = false;
  for (std::size_t reverse = bottom + 1; reverse-- > 0;) {
    const std::size_t level = reverse;
    if (buoyancy[level] > 1e-12) found_lfc = true;
    if (found_lfc && level > 0 && buoyancy[level] > 0.0 && buoyancy[level - 1] <= 0.0) {
      const Real fraction = buoyancy[level] / (buoyancy[level] - buoyancy[level - 1]);
      participation[level - 1] = fraction;
      for (std::size_t upper = 0; upper + 1 < level; ++upper)
        participation[upper] = 0.0;
      diagnostics.convection_top_pressure_pa =
          std::exp(std::log(input.pressure_full_pa[level]) +
                   fraction * std::log(input.pressure_full_pa[level - 1] /
                                       input.pressure_full_pa[level]));
      diagnostics.reaches_model_top = false;
      break;
    }
  }
}

}  // namespace

void simple_betts_miller_adjustment(const SimpleBettsMillerInput& input,
                                    SimpleBettsMillerResult& result) {
  input.thermodynamics.validate();
  const std::size_t levels = input.temperature_k.size();
  if (levels < 2 || input.vapor_mixing_ratio.size() != levels ||
      input.air_mass_kg_m2.size() != levels ||
      input.pressure_full_pa.size() != levels ||
      input.pressure_half_pa.size() != levels + 1)
    throw std::invalid_argument("Simple Betts-Miller column shape mismatch");
  require_positive(input.gravity_m_s2, "gravity");
  require_positive(input.relative_humidity_reference, "reference relative humidity");
  if (input.relative_humidity_reference > 1.0)
    throw std::invalid_argument("reference relative humidity must not exceed one");
  require_positive(input.relaxation_time_s, "convection relaxation time");
  require_positive(input.time_step_s, "convection time step");
  require_positive(input.minimum_temperature_k, "minimum convection temperature");
  for (std::size_t level = 0; level < levels; ++level) {
    validate_dilute_moist_state(input.temperature_k[level],
                                input.pressure_full_pa[level],
                                input.vapor_mixing_ratio[level], input.thermodynamics);
    require_positive(input.air_mass_kg_m2[level], "air mass");
    require_positive(input.pressure_half_pa[level], "half-level pressure");
    if (level > 0 &&
        !(input.pressure_full_pa[level] > input.pressure_full_pa[level - 1]))
      throw std::invalid_argument("full-level pressure must increase downward");
    if (!(input.pressure_full_pa[level] > input.pressure_half_pa[level]) ||
        !(input.pressure_full_pa[level] < input.pressure_half_pa[level + 1]))
      throw std::invalid_argument("full-level pressure must lie within its layer");
  }
  require_positive(input.pressure_half_pa.back(), "half-level pressure");

  result.temperature_k.assign(input.temperature_k.begin(), input.temperature_k.end());
  result.vapor_mixing_ratio.assign(input.vapor_mixing_ratio.begin(),
                                   input.vapor_mixing_ratio.end());
  result.reference_temperature_k.assign(input.temperature_k.begin(),
                                        input.temperature_k.end());
  result.reference_vapor_mixing_ratio.assign(input.vapor_mixing_ratio.begin(),
                                             input.vapor_mixing_ratio.end());
  result.diagnostics = {};
  const auto parcel = lift_surface_parcel(input);
  result.parcel_temperature_k = parcel.temperature;
  result.diagnostics.lcl_pressure_pa = parcel.lcl_pressure_pa;
  diagnose_buoyancy(input, result.parcel_temperature_k, result.participation_fraction,
                    result.diagnostics);
  if (!(result.diagnostics.cape_j_kg > 1e-12)) {
    result.diagnostics.reason = MoistConvectionReason::kNoCape;
    return;
  }

  Real water_excess = 0.0;
  Real thermal_sum = 0.0;
  Real effective_mass = 0.0;
  for (std::size_t level = 0; level < levels; ++level) {
    if (result.participation_fraction[level] == 0.0) continue;
    result.reference_temperature_k[level] = result.parcel_temperature_k[level];
    result.reference_vapor_mixing_ratio[level] =
        input.relative_humidity_reference *
        saturation_mixing_ratio(result.parcel_temperature_k[level],
                                input.pressure_full_pa[level], input.thermodynamics);
    const Real effective_layer_mass =
        input.air_mass_kg_m2[level] * result.participation_fraction[level];
    effective_mass += effective_layer_mass;
    water_excess += effective_layer_mass * (input.vapor_mixing_ratio[level] -
                                            result.reference_vapor_mixing_ratio[level]);
    thermal_sum += effective_layer_mass *
                   (result.reference_temperature_k[level] - input.temperature_k[level]);
  }
  const Real thermal_deficit = input.thermodynamics.heat_capacity_cp_j_kg_k /
                               input.thermodynamics.latent_heat_vaporization_j_kg *
                               thermal_sum;
  result.diagnostics.raw_water_excess_kg_m2 = water_excess;
  result.diagnostics.raw_thermal_deficit_kg_m2 = thermal_deficit;
  if (!(thermal_deficit > 0.0)) {
    result.diagnostics.reason = MoistConvectionReason::kThermalCondition;
    return;
  }

  if (water_excess > 0.0) {
    result.diagnostics.branch = MoistConvectionBranch::kDeep;
    const Real offset = input.thermodynamics.latent_heat_vaporization_j_kg *
                        (water_excess - thermal_deficit) /
                        (input.thermodynamics.heat_capacity_cp_j_kg_k * effective_mass);
    for (std::size_t level = 0; level < levels; ++level)
      if (result.participation_fraction[level] > 0.0)
        result.reference_temperature_k[level] += offset;
  } else {
    Real cumulative = 0.0;
    bool found = false;
    std::vector<Real> shallow_fraction(levels, 0.0);
    for (std::size_t reverse = levels; reverse > 0; --reverse) {
      const std::size_t level = reverse - 1;
      if (result.participation_fraction[level] == 0.0) continue;
      const Real full_layer_contribution =
          input.air_mass_kg_m2[level] * (result.reference_vapor_mixing_ratio[level] -
                                         input.vapor_mixing_ratio[level]);
      const Real maximum_fraction = result.participation_fraction[level];
      const Real next = cumulative + maximum_fraction * full_layer_contribution;
      if (cumulative != 0.0 && full_layer_contribution != 0.0 &&
          ((cumulative < 0.0 && next >= 0.0) || (cumulative > 0.0 && next <= 0.0))) {
        shallow_fraction[level] =
            std::clamp(-cumulative / full_layer_contribution, 0.0, maximum_fraction);
        result.diagnostics.convection_top_pressure_pa = std::exp(
            std::log(input.pressure_half_pa[level + 1]) +
            shallow_fraction[level] * std::log(input.pressure_half_pa[level] /
                                               input.pressure_half_pa[level + 1]));
        result.diagnostics.reaches_model_top = false;
        found = true;
        break;
      }
      shallow_fraction[level] = maximum_fraction;
      cumulative = next;
    }
    if (!found) {
      result.diagnostics.reason = MoistConvectionReason::kNoShallowSolution;
      result.participation_fraction.assign(levels, 0.0);
      return;
    }
    result.participation_fraction = shallow_fraction;
    Real mass = 0.0;
    Real temperature_difference = 0.0;
    for (std::size_t level = 0; level < levels; ++level) {
      const Real effective_layer_mass =
          input.air_mass_kg_m2[level] * result.participation_fraction[level];
      mass += effective_layer_mass;
      temperature_difference +=
          effective_layer_mass *
          (input.temperature_k[level] - result.reference_temperature_k[level]);
    }
    const Real offset = temperature_difference / mass;
    for (std::size_t level = 0; level < levels; ++level) {
      if (result.participation_fraction[level] > 0.0) {
        result.reference_temperature_k[level] += offset;
      } else {
        result.reference_temperature_k[level] = input.temperature_k[level];
        result.reference_vapor_mixing_ratio[level] = input.vapor_mixing_ratio[level];
      }
    }
    result.diagnostics.branch = MoistConvectionBranch::kShallow;
  }

  const Real relaxation =
      input.time_step_s / (input.relaxation_time_s + input.time_step_s);
  Real water_change = 0.0;
  Real enthalpy_change = 0.0;
  for (std::size_t level = 0; level < levels; ++level) {
    const Real fraction = result.participation_fraction[level];
    const Real temperature_increment =
        relaxation * fraction *
        (result.reference_temperature_k[level] - input.temperature_k[level]);
    const Real vapor_increment =
        relaxation * fraction *
        (result.reference_vapor_mixing_ratio[level] - input.vapor_mixing_ratio[level]);
    result.temperature_k[level] += temperature_increment;
    result.vapor_mixing_ratio[level] += vapor_increment;
    if (result.temperature_k[level] < input.minimum_temperature_k ||
        result.vapor_mixing_ratio[level] < 0.0 ||
        !std::isfinite(result.temperature_k[level]) ||
        !std::isfinite(result.vapor_mixing_ratio[level]))
      throw std::runtime_error("Simple Betts-Miller adjustment is inadmissible");
    const Real layer_water_change = input.air_mass_kg_m2[level] * vapor_increment;
    water_change += layer_water_change;
    enthalpy_change +=
        input.air_mass_kg_m2[level] *
        (input.thermodynamics.heat_capacity_cp_j_kg_k * temperature_increment +
         input.thermodynamics.latent_heat_vaporization_j_kg * vapor_increment);
    result.diagnostics.maximum_temperature_increment_k =
        std::max(result.diagnostics.maximum_temperature_increment_k,
                 std::abs(temperature_increment));
  }
  result.diagnostics.reason = MoistConvectionReason::kAdjusted;
  result.diagnostics.column_water_change_kg_m2 = water_change;
  result.diagnostics.convective_rain_kg_m2 =
      result.diagnostics.branch == MoistConvectionBranch::kDeep
          ? std::max(0.0, -water_change)
          : 0.0;
  result.diagnostics.moist_enthalpy_change_j_m2 = enthalpy_change;
}

void diagnose_sbm_reference(const SimpleBettsMillerInput& input,
                            SimpleBettsMillerReference& reference) {
  // Reuse the Phase 13 reference construction while the public split is introduced.
  // The finite relaxation is O(K) and is discarded here; the expensive parcel ascent
  // and buoyancy diagnosis are represented only by this call boundary.
  SimpleBettsMillerResult diagnosed;
  simple_betts_miller_adjustment(input, diagnosed);
  reference.parcel_temperature_k = std::move(diagnosed.parcel_temperature_k);
  reference.reference_temperature_k = std::move(diagnosed.reference_temperature_k);
  reference.reference_vapor_mixing_ratio =
      std::move(diagnosed.reference_vapor_mixing_ratio);
  reference.participation_fraction = std::move(diagnosed.participation_fraction);
  reference.diagnostics = diagnosed.diagnostics;
  reference.diagnostics.convective_rain_kg_m2 = 0.0;
  reference.diagnostics.column_water_change_kg_m2 = 0.0;
  reference.diagnostics.moist_enthalpy_change_j_m2 = 0.0;
  reference.diagnostics.maximum_temperature_increment_k = 0.0;
}

void apply_sbm_relaxation(const SimpleBettsMillerInput& input,
                          const SimpleBettsMillerReference& reference,
                          SimpleBettsMillerResult& result) {
  input.thermodynamics.validate();
  const std::size_t levels = input.temperature_k.size();
  if (levels < 2 || input.vapor_mixing_ratio.size() != levels ||
      input.air_mass_kg_m2.size() != levels ||
      reference.reference_temperature_k.size() != levels ||
      reference.reference_vapor_mixing_ratio.size() != levels ||
      reference.parcel_temperature_k.size() != levels ||
      reference.participation_fraction.size() != levels)
    throw std::invalid_argument("Simple Betts-Miller reference shape mismatch");
  require_positive(input.relaxation_time_s, "convection relaxation time");
  require_positive(input.time_step_s, "convection time step");
  require_positive(input.minimum_temperature_k, "minimum convection temperature");

  result.temperature_k.assign(input.temperature_k.begin(), input.temperature_k.end());
  result.vapor_mixing_ratio.assign(input.vapor_mixing_ratio.begin(),
                                   input.vapor_mixing_ratio.end());
  result.parcel_temperature_k = reference.parcel_temperature_k;
  result.reference_temperature_k = reference.reference_temperature_k;
  result.reference_vapor_mixing_ratio = reference.reference_vapor_mixing_ratio;
  result.participation_fraction = reference.participation_fraction;
  result.diagnostics = reference.diagnostics;
  if (reference.diagnostics.reason != MoistConvectionReason::kAdjusted) return;

  const Real relaxation =
      input.time_step_s / (input.relaxation_time_s + input.time_step_s);
  Real water_change = 0.0;
  Real enthalpy_change = 0.0;
  for (std::size_t level = 0; level < levels; ++level) {
    validate_dilute_moist_state(input.temperature_k[level],
                                input.pressure_full_pa[level],
                                input.vapor_mixing_ratio[level], input.thermodynamics);
    require_positive(input.air_mass_kg_m2[level], "air mass");
    const Real fraction = reference.participation_fraction[level];
    const Real temperature_increment =
        relaxation * fraction *
        (reference.reference_temperature_k[level] - input.temperature_k[level]);
    const Real vapor_increment = relaxation * fraction *
                                 (reference.reference_vapor_mixing_ratio[level] -
                                  input.vapor_mixing_ratio[level]);
    result.temperature_k[level] += temperature_increment;
    result.vapor_mixing_ratio[level] += vapor_increment;
    if (result.temperature_k[level] < input.minimum_temperature_k ||
        result.vapor_mixing_ratio[level] < 0.0 ||
        !std::isfinite(result.temperature_k[level]) ||
        !std::isfinite(result.vapor_mixing_ratio[level]))
      throw std::runtime_error("Simple Betts-Miller adjustment is inadmissible");
    const Real layer_water_change = input.air_mass_kg_m2[level] * vapor_increment;
    water_change += layer_water_change;
    enthalpy_change +=
        input.air_mass_kg_m2[level] *
        (input.thermodynamics.heat_capacity_cp_j_kg_k * temperature_increment +
         input.thermodynamics.latent_heat_vaporization_j_kg * vapor_increment);
    result.diagnostics.maximum_temperature_increment_k =
        std::max(result.diagnostics.maximum_temperature_increment_k,
                 std::abs(temperature_increment));
  }
  result.diagnostics.column_water_change_kg_m2 = water_change;
  result.diagnostics.convective_rain_kg_m2 =
      result.diagnostics.branch == MoistConvectionBranch::kDeep
          ? std::max(0.0, -water_change)
          : 0.0;
  result.diagnostics.moist_enthalpy_change_j_m2 = enthalpy_change;
}

SimpleBettsMillerResult simple_betts_miller_adjustment(
    const SimpleBettsMillerInput& input) {
  SimpleBettsMillerResult result;
  simple_betts_miller_adjustment(input, result);
  return result;
}

}  // namespace mps
