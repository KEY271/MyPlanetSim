#pragma once

#include <span>
#include <vector>

#include "myplanetsim/config/experiment_config.hpp"
#include "myplanetsim/core/planet_parameters.hpp"
#include "myplanetsim/grid/cubed_sphere_grid.hpp"
#include "myplanetsim/vertical/hybrid_pressure_coordinate.hpp"
#include "myplanetsim/vertical/hydrostatic_column.hpp"

namespace mps {

struct DryHydrostaticReferenceColumn {
  Real surface_pressure_pa = 0.0;
  Real temperature_k = 0.0;
  HybridPressureGeometry geometry;
  std::vector<Real> potential_temperature_k;
  std::vector<Real> potential_temperature_mass_k_kg_m2;
  HydrostaticColumn hydrostatic;
};

// Eigenvectors are stored mode-major. P10.06 registers the external mode; P10.09
// extends this same representation with the internal gravity modes.
struct DryHydrostaticVerticalModes {
  std::size_t levels = 0;
  std::vector<Real> phase_speed_m_s;
  std::vector<Real> eigenvalue_m2_s2;
  std::vector<Real> eigenvectors;

  [[nodiscard]] std::size_t mode_count() const noexcept {
    return phase_speed_m_s.size();
  }
  [[nodiscard]] std::span<const Real> eigenvector(std::size_t mode) const;
};

[[nodiscard]] DryHydrostaticReferenceColumn make_dry_hydrostatic_reference_column(
    const AtmosphericHybridCoordinate& coordinate, const PlanetParameters& planet,
    const SemiImplicitParameters& parameters);

[[nodiscard]] DryHydrostaticVerticalModes make_dry_hydrostatic_external_mode(
    const DryHydrostaticReferenceColumn& reference, const PlanetParameters& planet);

[[nodiscard]] std::vector<std::size_t> select_implicit_vertical_modes(
    const CubedSphereGrid& grid, const DryHydrostaticVerticalModes& modes,
    Real time_step_s, Real wave_cfl_threshold, std::size_t maximum_implicit_modes);

void project_onto_vertical_modes(const DryHydrostaticVerticalModes& modes,
                                 std::span<const Real> level_values,
                                 std::span<Real> mode_values);
void reconstruct_from_vertical_modes(const DryHydrostaticVerticalModes& modes,
                                     std::span<const Real> mode_values,
                                     std::span<Real> level_values);

}  // namespace mps
