#include "myplanetsim/dynamics/dry_hydrostatic_diffusion.hpp"

#include <stdexcept>

#include "myplanetsim/core/validation.hpp"
#include "myplanetsim/numerics/spherical_operators.hpp"

namespace mps {

DryHydrostaticDiffusionTendency dry_hydrostatic_diffusion_tendency(
    const CubedSphereGrid& grid, const DryHydrostaticDerived& derived,
    const DiffusionKind kind, const Real diffusion_coefficient) {
  require_non_negative(diffusion_coefficient, "dry_hydrostatic.diffusion_coefficient");
  const std::size_t cells = derived.cells;
  const std::size_t levels = derived.levels;
  if (cells != grid.cell_count()) {
    throw std::invalid_argument("dry diffusion cell count does not match grid");
  }
  DryHydrostaticDiffusionTendency tendency{
      .momentum = std::vector<Vec3>(cells * levels),
      .potential_temperature_mass = std::vector<Real>(cells * levels),
      .tracer_mass = std::vector<Real>(cells * levels),
  };
  if (kind == DiffusionKind::kNone) {
    if (diffusion_coefficient != 0.0) {
      throw std::invalid_argument("disabled diffusion requires zero coefficient");
    }
    return tendency;
  }
  require_positive(diffusion_coefficient, "active diffusion coefficient");
  // Biharmonic is the Laplacian applied twice with the opposite sign, matching the
  // shallow-water convention so one configured coefficient means the same thing.
  const bool biharmonic = kind == DiffusionKind::kBiharmonic;
  const Real signed_coefficient =
      biharmonic ? -diffusion_coefficient : diffusion_coefficient;

  std::vector<Vec3> velocity(cells);
  std::vector<Real> potential_temperature(cells);
  std::vector<Real> tracer(cells);
  for (std::size_t level = 0; level < levels; ++level) {
    for (std::size_t cell = 0; cell < cells; ++cell) {
      const auto n = dry_hydrostatic_offset(cell, level, levels);
      velocity[cell] = derived.velocity_m_s[n];
      potential_temperature[cell] = derived.potential_temperature_k[n];
      tracer[cell] = derived.tracer_mixing_ratio[n];
    }
    auto velocity_operator = finite_volume_vector_laplacian(grid, velocity);
    auto temperature_operator = finite_volume_laplacian(grid, potential_temperature);
    auto tracer_operator = finite_volume_laplacian(grid, tracer);
    if (biharmonic) {
      velocity_operator = finite_volume_vector_laplacian(grid, velocity_operator);
      temperature_operator = finite_volume_laplacian(grid, temperature_operator);
      tracer_operator = finite_volume_laplacian(grid, tracer_operator);
    }
    for (std::size_t cell = 0; cell < cells; ++cell) {
      const auto n = dry_hydrostatic_offset(cell, level, levels);
      const Real mass = derived.air_mass_kg_m2[n];
      const Vec3 momentum =
          signed_coefficient * mass *
          project_tangent(velocity_operator[cell], grid.cells()[cell].center);
      tendency.momentum[n] = momentum;
      tendency.potential_temperature_mass[n] =
          signed_coefficient * mass * temperature_operator[cell];
      tendency.tracer_mass[n] = signed_coefficient * mass * tracer_operator[cell];
      tendency.kinetic_energy_rate_w +=
          grid.cells()[cell].area_m2 * dot(velocity[cell], momentum);
    }
  }
  return tendency;
}

}  // namespace mps
