#include "myplanetsim/dynamics/dry_hydrostatic_fast_modes.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "myplanetsim/core/validation.hpp"

namespace mps {

std::span<const Real> DryHydrostaticVerticalModes::eigenvector(
    const std::size_t mode) const {
  if (mode >= mode_count() || eigenvectors.size() != mode_count() * levels)
    throw std::out_of_range("dry vertical mode index or shape is invalid");
  return {eigenvectors.data() + mode * levels, levels};
}

DryHydrostaticReferenceColumn make_dry_hydrostatic_reference_column(
    const AtmosphericHybridCoordinate& coordinate, const PlanetParameters& planet,
    const SemiImplicitParameters& parameters) {
  planet.validate();
  require_positive(parameters.reference_surface_pressure_pa,
                   "semi-implicit reference surface pressure");
  require_positive(parameters.reference_temperature_k,
                   "semi-implicit reference temperature");
  DryHydrostaticReferenceColumn result{
      .surface_pressure_pa = parameters.reference_surface_pressure_pa,
      .temperature_k = parameters.reference_temperature_k,
      .geometry = coordinate.geometry(parameters.reference_surface_pressure_pa,
                                      planet.gravity_m_s2, planet.gas_constant_j_kg_k,
                                      planet.heat_capacity_cp_j_kg_k,
                                      planet.reference_pressure_pa)};
  const auto levels = coordinate.levels();
  result.potential_temperature_k.resize(levels);
  result.potential_temperature_mass_k_kg_m2.resize(levels);
  for (std::size_t level = 0; level < levels; ++level) {
    const Real potential_temperature =
        parameters.reference_temperature_k / result.geometry.exner_full[level];
    result.potential_temperature_k[level] = potential_temperature;
    result.potential_temperature_mass_k_kg_m2[level] =
        result.geometry.air_mass_kg_m2[level] * potential_temperature;
  }
  result.hydrostatic = integrate_hydrostatic_column(
      result.geometry, result.potential_temperature_k, planet.heat_capacity_cp_j_kg_k,
      planet.gravity_m_s2, 0.0);
  return result;
}

DryHydrostaticVerticalModes make_dry_hydrostatic_external_mode(
    const DryHydrostaticReferenceColumn& reference, const PlanetParameters& planet) {
  planet.validate();
  const std::size_t levels = reference.geometry.air_mass_kg_m2.size();
  if (levels == 0 || reference.potential_temperature_k.size() != levels ||
      reference.temperature_k <= 0.0 || !std::isfinite(reference.temperature_k))
    throw std::invalid_argument("dry reference column shape or temperature is invalid");
  const Real gamma = planet.heat_capacity_cp_j_kg_k /
                     (planet.heat_capacity_cp_j_kg_k - planet.gas_constant_j_kg_k);
  const Real eigenvalue = gamma * planet.gas_constant_j_kg_k * reference.temperature_k;
  const Real normalization = 1.0 / std::sqrt(static_cast<Real>(levels));
  return {.levels = levels,
          .phase_speed_m_s = {std::sqrt(eigenvalue)},
          .eigenvalue_m2_s2 = {eigenvalue},
          .eigenvectors = std::vector<Real>(levels, normalization)};
}

std::vector<std::size_t> select_implicit_vertical_modes(
    const CubedSphereGrid& grid, const DryHydrostaticVerticalModes& modes,
    const Real time_step_s, const Real wave_cfl_threshold,
    const std::size_t maximum_implicit_modes) {
  require_positive(time_step_s, "semi-implicit mode-selection time step");
  require_positive(wave_cfl_threshold, "semi-implicit wave CFL threshold");
  if (maximum_implicit_modes == 0)
    throw std::invalid_argument("maximum implicit mode count must be positive");
  if (modes.mode_count() == 0 || modes.eigenvalue_m2_s2.size() != modes.mode_count() ||
      modes.eigenvectors.size() != modes.mode_count() * modes.levels)
    throw std::invalid_argument("dry vertical modes are empty or malformed");
  Real maximum_geometry_factor_m_inverse = 0.0;
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
    Real perimeter = 0.0;
    for (const auto& edge : grid.cell_cache()[cell].edges)
      perimeter += grid.edges()[edge.edge].length_m;
    maximum_geometry_factor_m_inverse = std::max(
        maximum_geometry_factor_m_inverse, perimeter / grid.cells()[cell].area_m2);
  }
  std::vector<std::size_t> selected;
  Real previous_speed = std::numeric_limits<Real>::infinity();
  for (std::size_t mode = 0; mode < modes.mode_count(); ++mode) {
    const Real speed = modes.phase_speed_m_s[mode];
    if (!(speed > 0.0) || !std::isfinite(speed) || speed > previous_speed)
      throw std::invalid_argument(
          "dry vertical mode speeds must be positive and descending");
    previous_speed = speed;
    const Real courant = time_step_s * speed * maximum_geometry_factor_m_inverse;
    if (courant > wave_cfl_threshold) selected.push_back(mode);
  }
  if (selected.size() > maximum_implicit_modes)
    throw std::runtime_error(
        "required implicit vertical modes exceed the configured maximum");
  return selected;
}

void project_onto_vertical_modes(const DryHydrostaticVerticalModes& modes,
                                 const std::span<const Real> level_values,
                                 const std::span<Real> mode_values) {
  if (level_values.size() != modes.levels || mode_values.size() != modes.mode_count())
    throw std::invalid_argument("vertical mode projection shapes differ");
  for (std::size_t mode = 0; mode < modes.mode_count(); ++mode) {
    mode_values[mode] = 0.0;
    const auto eigenvector = modes.eigenvector(mode);
    for (std::size_t level = 0; level < modes.levels; ++level)
      mode_values[mode] += eigenvector[level] * level_values[level];
  }
}

void reconstruct_from_vertical_modes(const DryHydrostaticVerticalModes& modes,
                                     const std::span<const Real> mode_values,
                                     const std::span<Real> level_values) {
  if (level_values.size() != modes.levels || mode_values.size() != modes.mode_count())
    throw std::invalid_argument("vertical mode reconstruction shapes differ");
  std::fill(level_values.begin(), level_values.end(), 0.0);
  for (std::size_t mode = 0; mode < modes.mode_count(); ++mode) {
    const auto eigenvector = modes.eigenvector(mode);
    for (std::size_t level = 0; level < modes.levels; ++level)
      level_values[level] += mode_values[mode] * eigenvector[level];
  }
}

}  // namespace mps
