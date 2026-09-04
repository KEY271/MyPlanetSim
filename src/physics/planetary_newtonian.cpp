#include "myplanetsim/physics/planetary_newtonian.hpp"

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

}  // namespace

HeldSuarezRates planetary_newtonian_rates(const Real latitude_rad, const Real chi,
                                          const Real pressure_pa,
                                          const Real surface_pressure_pa,
                                          const PlanetParameters& planet) {
  if (!std::isfinite(latitude_rad) || !std::isfinite(chi) || chi < -1.0 || chi > 1.0 ||
      !(pressure_pa > 0.0) || !(surface_pressure_pa > 0.0) ||
      pressure_pa > surface_pressure_pa)
    throw std::invalid_argument("planetary Newtonian forcing input is invalid");
  const Real cosine = std::cos(latitude_rad);
  const Real pressure_ratio = pressure_pa / planet.reference_pressure_pa;
  const Real equilibrium = std::max(
      200.0, (315.0 + 60.0 * chi - 10.0 * std::log(pressure_ratio) * cosine * cosine) *
                 std::pow(pressure_ratio, planet.kappa()));
  const Real boundary_weight =
      std::max(0.0, (pressure_pa / surface_pressure_pa - 0.7) / 0.3);
  const Real atmospheric_rate = 1.0 / (40.0 * kSecondsPerDay);
  const Real surface_rate = 1.0 / (4.0 * kSecondsPerDay);
  return {.equilibrium_temperature_k = equilibrium,
          .temperature_relaxation_rate_s_1 =
              atmospheric_rate +
              (surface_rate - atmospheric_rate) * boundary_weight * std::pow(cosine, 4),
          .rayleigh_drag_rate_s_1 = boundary_weight / kSecondsPerDay};
}

HeldSuarezTendency planetary_newtonian_tendency(
    const CubedSphereGrid& grid, const AtmosphericHybridCoordinate& coordinate,
    const DryHydrostaticDerived& derived,
    const std::span<const Real> surface_pressure_pa, const PlanetParameters& planet,
    const ForcingGeometry geometry, const OrbitState* orbit_state) {
  if (geometry == ForcingGeometry::kAxisymmetric)
    return held_suarez_tendency(grid, coordinate, derived, surface_pressure_pa, planet);
  if (orbit_state == nullptr)
    throw std::invalid_argument("substellar forcing requires an orbit state");
  const std::size_t volume = derived.cells * derived.levels;
  if (derived.cells != grid.cell_count() || derived.levels == 0 ||
      derived.pressure_pa.size() != volume || derived.air_mass_kg_m2.size() != volume ||
      derived.velocity_m_s.size() != volume || derived.temperature_k.size() != volume ||
      surface_pressure_pa.size() != derived.cells ||
      coordinate.levels() != derived.levels)
    throw std::invalid_argument("planetary Newtonian derived-state shape mismatch");

  HeldSuarezTendency result;
  result.potential_temperature_mass_k_kg_m2_s.resize(volume);
  result.horizontal_momentum_mass_kg_m_s2.resize(volume);
  result.equilibrium_temperature_k.resize(volume);
  result.temperature_relaxation_rate_s_1.resize(volume);
  result.rayleigh_drag_rate_s_1.resize(volume);
  const Real cv = planet.heat_capacity_cv_j_kg_k();
  for (std::size_t cell = 0; cell < derived.cells; ++cell) {
    const Vec3 position = grid.cells()[cell].center;
    const Real latitude = std::asin(position.z);
    const Real chi = dot(position, orbit_state->star_direction_planet_fixed);
    const Vec3 east = eastward_unit(position);
    const Vec3 north = cross(position, east);
    const Real area = grid.cells()[cell].area_m2;
    const auto vertical_geometry = coordinate.geometry(
        surface_pressure_pa[cell], planet.gravity_m_s2, planet.gas_constant_j_kg_k,
        planet.heat_capacity_cp_j_kg_k, planet.reference_pressure_pa);
    std::vector<Real> theta_mass_rates(derived.levels);
    for (std::size_t level = 0; level < derived.levels; ++level) {
      const auto offset = dry_hydrostatic_offset(cell, level, derived.levels);
      const auto rates =
          planetary_newtonian_rates(latitude, chi, derived.pressure_pa[offset],
                                    surface_pressure_pa[cell], planet);
      const Real temperature_rate =
          -rates.temperature_relaxation_rate_s_1 *
          (derived.temperature_k[offset] - rates.equilibrium_temperature_k);
      const Real exner = std::pow(
          derived.pressure_pa[offset] / planet.reference_pressure_pa, planet.kappa());
      const Real theta_mass_rate =
          derived.air_mass_kg_m2[offset] * temperature_rate / exner;
      const Vec3 momentum_rate = -derived.air_mass_kg_m2[offset] *
                                 rates.rayleigh_drag_rate_s_1 *
                                 derived.velocity_m_s[offset];
      result.potential_temperature_mass_k_kg_m2_s[offset] = theta_mass_rate;
      theta_mass_rates[level] = theta_mass_rate;
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
    for (std::size_t reverse = derived.levels; reverse > 0; --reverse) {
      const std::size_t level = reverse - 1;
      const auto offset = dry_hydrostatic_offset(cell, level, derived.levels);
      const Real theta_rate = theta_mass_rates[level] / derived.air_mass_kg_m2[offset];
      const Real geopotential_full_rate =
          geopotential_half_rate + planet.heat_capacity_cp_j_kg_k * theta_rate *
                                       (vertical_geometry.exner_half[level + 1] -
                                        vertical_geometry.exner_full[level]);
      result.diagnostics.thermal_energy_rate_w +=
          area * derived.air_mass_kg_m2[offset] * geopotential_full_rate;
      geopotential_half_rate += planet.heat_capacity_cp_j_kg_k * theta_rate *
                                (vertical_geometry.exner_half[level + 1] -
                                 vertical_geometry.exner_half[level]);
    }
  }
  return result;
}

}  // namespace mps
