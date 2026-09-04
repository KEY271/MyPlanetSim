#include "myplanetsim/dynamics/dry_hydrostatic_reconstruction.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "myplanetsim/numerics/spherical_operators.hpp"

namespace mps {
namespace {

struct ScalarReconstruction {
  std::vector<Vec3> gradient;
  std::vector<Real> factor;
};

[[nodiscard]] Vec3 logarithmic_displacement(const Vec3 from, const Vec3 to,
                                            const Real radius_m) {
  const Real cosine = std::clamp(dot(from, to), -1.0, 1.0);
  const Real angle = std::acos(cosine);
  const Real sine = std::sin(angle);
  if (!(sine > 0.0)) return {};
  return (radius_m * angle / sine) * (to - cosine * from);
}

[[nodiscard]] Real barth_factor(const Real center, const Real increment,
                                const Real minimum, const Real maximum) {
  if (increment > 0.0) return std::min(1.0, (maximum - center) / increment);
  if (increment < 0.0) return std::min(1.0, (minimum - center) / increment);
  return 1.0;
}

[[nodiscard]] ScalarReconstruction prepare_scalar(const CubedSphereGrid& grid,
                                                  const std::vector<Real>& values,
                                                  const LimiterKind limiter) {
  ScalarReconstruction result{.gradient = least_squares_gradient(grid, values),
                              .factor = std::vector<Real>(values.size(), 1.0)};
  if (limiter == LimiterKind::kNone) return result;
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
    Real minimum = values[cell];
    Real maximum = values[cell];
    for (const auto edge_id : grid.cell_edges(grid.cell_id(cell))) {
      const auto neighbor =
          grid.cell_index(grid.neighbor_across(edge_id, grid.cell_id(cell)));
      minimum = std::min(minimum, values[neighbor]);
      maximum = std::max(maximum, values[neighbor]);
    }
    for (const auto edge_id : grid.cell_edges(grid.cell_id(cell))) {
      const auto displacement = logarithmic_displacement(
          grid.cells()[cell].center, grid.edge(edge_id).center, grid.radius_m());
      result.factor[cell] =
          std::min(result.factor[cell],
                   barth_factor(values[cell], dot(result.gradient[cell], displacement),
                                minimum, maximum));
    }
    result.factor[cell] = std::clamp(result.factor[cell], 0.0, 1.0);
  }
  return result;
}

[[nodiscard]] Real reconstruct_scalar(const CubedSphereGrid& grid,
                                      const std::size_t cell, const EdgeGeometry& edge,
                                      const std::vector<Real>& values,
                                      const ScalarReconstruction& prepared) {
  const auto displacement =
      logarithmic_displacement(grid.cells()[cell].center, edge.center, grid.radius_m());
  return values[cell] +
         prepared.factor[cell] * dot(prepared.gradient[cell], displacement);
}

[[nodiscard]] std::vector<Real> level_scalar(const std::vector<Real>& volume,
                                             const std::size_t cells,
                                             const std::size_t levels,
                                             const std::size_t level) {
  std::vector<Real> result(cells);
  for (std::size_t cell = 0; cell < cells; ++cell)
    result[cell] = volume[dry_hydrostatic_offset(cell, level, levels)];
  return result;
}

[[nodiscard]] std::vector<Vec3> level_vector(const std::vector<Vec3>& volume,
                                             const std::size_t cells,
                                             const std::size_t levels,
                                             const std::size_t level) {
  std::vector<Vec3> result(cells);
  for (std::size_t cell = 0; cell < cells; ++cell)
    result[cell] = volume[dry_hydrostatic_offset(cell, level, levels)];
  return result;
}

}  // namespace

const DryHydrostaticFaceStates& DryHydrostaticReconstruction::at(
    const std::size_t edge, const std::size_t level) const {
  if (level >= levels || edge * levels + level >= edge_levels.size())
    throw std::out_of_range("dry reconstruction edge-level is out of range");
  return edge_levels[edge * levels + level];
}

DryHydrostaticReconstruction reconstruct_dry_hydrostatic_face_states(
    const CubedSphereGrid& grid, const DryHydrostaticDerived& derived,
    const ReconstructionKind reconstruction, const LimiterKind limiter) {
  const auto cells = grid.cell_count();
  const auto levels = derived.levels;
  const auto volume = cells * levels;
  if (derived.cells != cells || levels == 0 ||
      derived.air_mass_kg_m2.size() != volume ||
      derived.velocity_m_s.size() != volume ||
      derived.potential_temperature_k.size() != volume ||
      derived.tracer_mixing_ratio.size() != volume ||
      derived.temperature_k.size() != volume)
    throw std::invalid_argument("dry reconstruction shape mismatch");

  DryHydrostaticReconstruction result{
      .levels = levels,
      .edge_levels = std::vector<DryHydrostaticFaceStates>(grid.edge_count() * levels),
      .limiter_activations = 0};
  for (std::size_t level = 0; level < levels; ++level) {
    const auto mass = level_scalar(derived.air_mass_kg_m2, cells, levels, level);
    const auto velocity = level_vector(derived.velocity_m_s, cells, levels, level);
    const auto theta =
        level_scalar(derived.potential_temperature_k, cells, levels, level);
    const auto tracer = level_scalar(derived.tracer_mixing_ratio, cells, levels, level);
    const auto temperature = level_scalar(derived.temperature_k, cells, levels, level);

    if (reconstruction == ReconstructionKind::kPiecewiseConstant) {
      for (const auto& edge : grid.edges()) {
        const auto left = grid.cell_index(edge.left_cell);
        const auto right = grid.cell_index(edge.right_cell);
        result.edge_levels[edge.id * levels + level] = {
            .left = {.air_mass_kg_m2 = mass[left],
                     .velocity_m_s = project_tangent(velocity[left], edge.center),
                     .potential_temperature_k = theta[left],
                     .tracer_mixing_ratio = tracer[left],
                     .temperature_k = temperature[left]},
            .right = {.air_mass_kg_m2 = mass[right],
                      .velocity_m_s = project_tangent(velocity[right], edge.center),
                      .potential_temperature_k = theta[right],
                      .tracer_mixing_ratio = tracer[right],
                      .temperature_k = temperature[right]}};
      }
      continue;
    }

    const auto mass_reconstruction = prepare_scalar(grid, mass, limiter);
    const auto theta_reconstruction = prepare_scalar(grid, theta, limiter);
    const auto tracer_reconstruction = prepare_scalar(grid, tracer, limiter);
    const auto temperature_reconstruction = prepare_scalar(grid, temperature, limiter);
    const auto velocity_gradient = least_squares_vector_gradient(grid, velocity);
    for (std::size_t cell = 0; cell < cells; ++cell) {
      if (mass_reconstruction.factor[cell] < 1.0 - 1.0e-14 ||
          theta_reconstruction.factor[cell] < 1.0 - 1.0e-14 ||
          tracer_reconstruction.factor[cell] < 1.0 - 1.0e-14 ||
          temperature_reconstruction.factor[cell] < 1.0 - 1.0e-14)
        ++result.limiter_activations;
    }
    for (const auto& edge : grid.edges()) {
      const auto left = grid.cell_index(edge.left_cell);
      const auto right = grid.cell_index(edge.right_cell);
      const auto face = [&](const std::size_t cell) {
        return DryHydrostaticPrimitive{
            .air_mass_kg_m2 =
                reconstruct_scalar(grid, cell, edge, mass, mass_reconstruction),
            .velocity_m_s = reconstruct_tangent_vector(
                grid, cell, velocity[cell], edge.center, velocity_gradient[cell]),
            .potential_temperature_k =
                reconstruct_scalar(grid, cell, edge, theta, theta_reconstruction),
            .tracer_mixing_ratio =
                reconstruct_scalar(grid, cell, edge, tracer, tracer_reconstruction),
            .temperature_k = reconstruct_scalar(grid, cell, edge, temperature,
                                                temperature_reconstruction)};
      };
      result.edge_levels[edge.id * levels + level] = {.left = face(left),
                                                      .right = face(right)};
    }
  }
  return result;
}

}  // namespace mps
