#include "myplanetsim/physics/held_suarez.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace mps {
namespace {

constexpr Real kSecondsPerDay = 86400.0;

[[nodiscard]] Vec3 eastward_unit(const Vec3 position) {
  const Vec3 east{-position.y, position.x, 0.0};
  const Real magnitude = norm(east);
  return magnitude > 0.0 ? east / magnitude : Vec3{0.0, 1.0, 0.0};
}

[[nodiscard]] HeldSuarezRates held_suarez_rates_from_trig(
    const Real sine, const Real cosine, const Real pressure_pa,
    const Real surface_pressure_pa, const PlanetParameters& planet,
    const Real exner = -1.0) {
  if (!std::isfinite(sine) || !std::isfinite(cosine) || !(pressure_pa > 0.0) ||
      !(surface_pressure_pa > 0.0) || pressure_pa > surface_pressure_pa) {
    throw std::invalid_argument("Held-Suarez pressure or latitude is invalid");
  }
  const Real pressure_ratio = pressure_pa / planet.reference_pressure_pa;
  const Real pressure_exner =
      exner > 0.0 ? exner : std::pow(pressure_ratio, planet.kappa());
  const Real cosine_squared = cosine * cosine;
  const Real equilibrium = std::max(
      200.0,
      (315.0 - 60.0 * sine * sine - 10.0 * std::log(pressure_ratio) * cosine_squared) *
          pressure_exner);
  const Real boundary_weight =
      std::max(0.0, (pressure_pa / surface_pressure_pa - 0.7) / 0.3);
  const Real atmospheric_rate = 1.0 / (40.0 * kSecondsPerDay);
  const Real surface_rate = 1.0 / (4.0 * kSecondsPerDay);
  return {.equilibrium_temperature_k = equilibrium,
          .temperature_relaxation_rate_s_1 =
              atmospheric_rate + (surface_rate - atmospheric_rate) * boundary_weight *
                                     cosine_squared * cosine_squared,
          .rayleigh_drag_rate_s_1 = boundary_weight / kSecondsPerDay};
}

}  // namespace

HeldSuarezRates held_suarez_rates(const Real latitude_rad, const Real pressure_pa,
                                  const Real surface_pressure_pa,
                                  const PlanetParameters& planet) {
  if (!std::isfinite(latitude_rad) || !(pressure_pa > 0.0) ||
      !(surface_pressure_pa > 0.0) || pressure_pa > surface_pressure_pa) {
    throw std::invalid_argument("Held-Suarez pressure or latitude is invalid");
  }
  return held_suarez_rates_from_trig(std::sin(latitude_rad), std::cos(latitude_rad),
                                     pressure_pa, surface_pressure_pa, planet);
}

HeldSuarezTendency held_suarez_tendency(const CubedSphereGrid& grid,
                                        const AtmosphericHybridCoordinate& coordinate,
                                        const DryHydrostaticDerived& derived,
                                        const std::span<const Real> surface_pressure_pa,
                                        const PlanetParameters& planet) {
  HeldSuarezTendency result;
  HeldSuarezWorkspace workspace;
  held_suarez_tendency(grid, coordinate, derived, surface_pressure_pa, planet, result,
                       workspace);
  return result;
}

void held_suarez_tendency(const CubedSphereGrid& grid,
                          const AtmosphericHybridCoordinate& coordinate,
                          const DryHydrostaticDerived& derived,
                          const std::span<const Real> surface_pressure_pa,
                          const PlanetParameters& planet, HeldSuarezTendency& result,
                          HeldSuarezWorkspace& workspace) {
  const std::size_t volume = derived.cells * derived.levels;
  if (derived.cells != grid.cell_count() || derived.levels == 0 ||
      derived.pressure_pa.size() != volume || derived.air_mass_kg_m2.size() != volume ||
      derived.velocity_m_s.size() != volume || derived.temperature_k.size() != volume) {
    throw std::invalid_argument("Held-Suarez derived-state shape mismatch");
  }
  if (surface_pressure_pa.size() != derived.cells) {
    throw std::invalid_argument("Held-Suarez surface-pressure shape mismatch");
  }
  if (coordinate.levels() != derived.levels) {
    throw std::invalid_argument("Held-Suarez vertical-coordinate shape mismatch");
  }

  result.diagnostics = {};
  result.potential_temperature_mass_k_kg_m2_s.resize(volume);
  result.horizontal_momentum_mass_kg_m_s2.resize(volume);
  result.equilibrium_temperature_k.resize(volume);
  result.temperature_relaxation_rate_s_1.resize(volume);
  result.rayleigh_drag_rate_s_1.resize(volume);
  workspace.theta_mass_rates.resize(derived.levels);
  const Real cv = planet.heat_capacity_cp_j_kg_k - planet.gas_constant_j_kg_k;
  const bool has_cached_exner =
      derived.exner_half.size() == derived.cells * (derived.levels + 1) &&
      derived.exner_full.size() == volume;

  for (std::size_t cell = 0; cell < derived.cells; ++cell) {
    const Vec3 position = grid.cells()[cell].center;
    const Real sine_latitude = position.z;
    const Real cosine_latitude =
        std::sqrt(std::max(0.0, 1.0 - sine_latitude * sine_latitude));
    const Vec3 east = eastward_unit(position);
    const Vec3 north = cross(position, east);
    const Real area = grid.cells()[cell].area_m2;
    if (!has_cached_exner) {
      coordinate.geometry(surface_pressure_pa[cell], planet.gravity_m_s2,
                          planet.gas_constant_j_kg_k, planet.heat_capacity_cp_j_kg_k,
                          planet.reference_pressure_pa, workspace.vertical_geometry);
    }
    const std::span<const Real> exner_half =
        has_cached_exner
            ? std::span<const Real>(
                  derived.exner_half.data() + cell * (derived.levels + 1),
                  derived.levels + 1)
            : std::span<const Real>(workspace.vertical_geometry.exner_half);
    const std::span<const Real> exner_full =
        has_cached_exner
            ? std::span<const Real>(derived.exner_full.data() + cell * derived.levels,
                                    derived.levels)
            : std::span<const Real>(workspace.vertical_geometry.exner_full);
    for (std::size_t level = 0; level < derived.levels; ++level) {
      const auto offset = dry_hydrostatic_offset(cell, level, derived.levels);
      const auto rates = held_suarez_rates_from_trig(
          sine_latitude, cosine_latitude, derived.pressure_pa[offset],
          surface_pressure_pa[cell], planet, derived.exner_full[offset]);
      const Real temperature_rate =
          -rates.temperature_relaxation_rate_s_1 *
          (derived.temperature_k[offset] - rates.equilibrium_temperature_k);
      const Real theta_mass_rate = derived.air_mass_kg_m2[offset] * temperature_rate /
                                   derived.exner_full[offset];
      const Vec3 momentum_rate = -derived.air_mass_kg_m2[offset] *
                                 rates.rayleigh_drag_rate_s_1 *
                                 derived.velocity_m_s[offset];

      result.potential_temperature_mass_k_kg_m2_s[offset] = theta_mass_rate;
      workspace.theta_mass_rates[level] = theta_mass_rate;
      result.horizontal_momentum_mass_kg_m_s2[offset] = momentum_rate;
      result.equilibrium_temperature_k[offset] = rates.equilibrium_temperature_k;
      result.temperature_relaxation_rate_s_1[offset] =
          rates.temperature_relaxation_rate_s_1;
      result.rayleigh_drag_rate_s_1[offset] = rates.rayleigh_drag_rate_s_1;
      result.diagnostics.potential_temperature_mass_rate_k_kg_s +=
          area * theta_mass_rate;
      result.diagnostics.eastward_momentum_rate_n += area * dot(east, momentum_rate);
      result.diagnostics.northward_momentum_rate_n += area * dot(north, momentum_rate);
      result.diagnostics.thermal_energy_rate_w +=
          area * cv * derived.air_mass_kg_m2[offset] * temperature_rate;
      result.diagnostics.rayleigh_drag_work_w +=
          area * dot(derived.velocity_m_s[offset], momentum_rate);
    }
    Real geopotential_half_rate = 0.0;
    Real column_potential_energy_rate = 0.0;
    for (std::size_t reverse = derived.levels; reverse > 0; --reverse) {
      const std::size_t level = reverse - 1;
      const auto offset = dry_hydrostatic_offset(cell, level, derived.levels);
      const Real theta_rate =
          workspace.theta_mass_rates[level] / derived.air_mass_kg_m2[offset];
      const Real geopotential_full_rate =
          geopotential_half_rate + planet.heat_capacity_cp_j_kg_k * theta_rate *
                                       (exner_half[level + 1] - exner_full[level]);
      column_potential_energy_rate +=
          derived.air_mass_kg_m2[offset] * geopotential_full_rate;
      geopotential_half_rate += planet.heat_capacity_cp_j_kg_k * theta_rate *
                                (exner_half[level + 1] - exner_half[level]);
    }
    result.diagnostics.thermal_energy_rate_w += area * column_potential_energy_rate;
  }
}

}  // namespace mps
