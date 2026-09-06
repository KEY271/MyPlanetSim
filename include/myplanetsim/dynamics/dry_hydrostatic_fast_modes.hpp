#pragma once

#include <iosfwd>
#include <span>
#include <vector>

#include "myplanetsim/config/experiment_config.hpp"
#include "myplanetsim/core/planet_parameters.hpp"
#include "myplanetsim/grid/cubed_sphere_grid.hpp"
#include "myplanetsim/vertical/hybrid_pressure_coordinate.hpp"
#include "myplanetsim/vertical/hydrostatic_column.hpp"

namespace mps {

struct DryHydrostaticFastOperator;

struct DryHydrostaticReferenceColumn {
  Real surface_pressure_pa = 0.0;
  Real temperature_k = 0.0;
  HybridPressureGeometry geometry;
  std::vector<Real> potential_temperature_k;
  std::vector<Real> potential_temperature_mass_k_kg_m2;
  HydrostaticColumn hydrostatic;
};

// Eigenvectors are stored mode-major. The full basis contains the external mode
// followed by internal gravity modes in descending phase-speed order.
struct DryHydrostaticVerticalModes {
  std::size_t levels = 0;
  std::vector<Real> phase_speed_m_s;
  std::vector<Real> eigenvalue_m2_s2;
  // Right eigenvectors are mode-major. inverse_eigenvectors contains the
  // corresponding rows of V^-1 and is used for projection when the vertical
  // structure matrix is not symmetric in prognostic momentum variables.
  std::vector<Real> eigenvectors;
  std::vector<Real> inverse_eigenvectors;
  std::vector<Real> vertical_structure_m2_s2;

  [[nodiscard]] std::size_t mode_count() const noexcept {
    return phase_speed_m_s.size();
  }
  [[nodiscard]] std::span<const Real> eigenvector(std::size_t mode) const;
};

[[nodiscard]] DryHydrostaticReferenceColumn make_dry_hydrostatic_reference_column(
    const AtmosphericHybridCoordinate& coordinate, const PlanetParameters& planet,
    const SemiImplicitParameters& parameters);

// Per-level form used by the per-step reference update. The column stays horizontally
// uniform, so the vertical structure matrix remains common to every cell and a single
// eigendecomposition per step is still sufficient. `temperature_k` is set to the
// mass-weighted mean of the profile and drives the external Lamb mode.
[[nodiscard]] DryHydrostaticReferenceColumn make_dry_hydrostatic_reference_column(
    const AtmosphericHybridCoordinate& coordinate, const PlanetParameters& planet,
    Real reference_surface_pressure_pa, std::span<const Real> temperature_profile_k);

[[nodiscard]] DryHydrostaticVerticalModes make_dry_hydrostatic_external_mode(
    const DryHydrostaticReferenceColumn& reference, const PlanetParameters& planet);

[[nodiscard]] DryHydrostaticVerticalModes make_dry_hydrostatic_vertical_modes(
    const DryHydrostaticReferenceColumn& reference, const PlanetParameters& planet,
    const DryHydrostaticFastOperator& fast_operator);

[[nodiscard]] std::vector<std::size_t> select_implicit_vertical_modes(
    const CubedSphereGrid& grid, const DryHydrostaticVerticalModes& modes,
    Real time_step_s, Real wave_cfl_threshold, std::size_t maximum_implicit_modes);

[[nodiscard]] Real maximum_vertical_mode_courant(
    const CubedSphereGrid& grid, const DryHydrostaticVerticalModes& modes,
    Real time_step_s);

void write_dry_hydrostatic_vertical_mode_metadata(
    std::ostream& output, const DryHydrostaticVerticalModes& modes);

void project_onto_vertical_modes(const DryHydrostaticVerticalModes& modes,
                                 std::span<const Real> level_values,
                                 std::span<Real> mode_values);
void reconstruct_from_vertical_modes(const DryHydrostaticVerticalModes& modes,
                                     std::span<const Real> mode_values,
                                     std::span<Real> level_values);

}  // namespace mps
