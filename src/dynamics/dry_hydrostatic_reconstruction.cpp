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

[[nodiscard]] Real barth_factor(const Real center, const Real increment,
                                const Real minimum, const Real maximum) {
  if (increment > 0.0) return std::min(1.0, (maximum - center) / increment);
  if (increment < 0.0) return std::min(1.0, (minimum - center) / increment);
  return 1.0;
}

// Shrinks a velocity increment until the reconstructed face speed stays inside the
// neighbourhood maximum, matching the shallow-water reconstruction.
[[nodiscard]] Real speed_factor(const Vec3 center, const Vec3 increment,
                                const Real maximum_speed) {
  if (norm(center + increment) <= maximum_speed) return 1.0;
  Real low = 0.0;
  Real high = 1.0;
  for (int iteration = 0; iteration < 54; ++iteration) {
    const Real middle = 0.5 * (low + high);
    if (norm(center + middle * increment) <= maximum_speed) {
      low = middle;
    } else {
      high = middle;
    }
  }
  return low;
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
    for (const auto& edge : grid.cell_cache()[cell].edges) {
      minimum = std::min(minimum, values[edge.neighbor]);
      maximum = std::max(maximum, values[edge.neighbor]);
    }
    for (const auto& edge : grid.cell_cache()[cell].edges) {
      result.factor[cell] =
          std::min(result.factor[cell],
                   barth_factor(values[cell],
                                dot(result.gradient[cell], edge.face_displacement_m),
                                minimum, maximum));
    }
    result.factor[cell] = std::clamp(result.factor[cell], 0.0, 1.0);
  }
  return result;
}

[[nodiscard]] Real reconstruct_scalar(const std::size_t cell,
                                      const CachedCellEdgeGeometry& edge,
                                      const std::vector<Real>& values,
                                      const ScalarReconstruction& prepared) {
  return values[cell] +
         prepared.factor[cell] * dot(prepared.gradient[cell], edge.face_displacement_m);
}

[[nodiscard]] const CachedCellEdgeGeometry& cached_cell_edge(
    const CubedSphereGrid& grid, const std::size_t cell, const std::size_t edge) {
  for (const auto& cached : grid.cell_cache()[cell].edges) {
    if (cached.edge == edge) return cached;
  }
  throw std::logic_error("cached edge is not incident to cell");
}

struct VelocityReconstruction {
  std::vector<TangentVectorGradient> gradient;
  std::vector<Real> factor;
};

// Limits the reconstructed tangent velocity the same way the shallow-water path does:
// the face normal and tangent components stay inside the neighbourhood range and the
// face speed stays inside the neighbourhood maximum. Without this the unlimited
// least-squares vector gradient overshoots at cubed-sphere seams.
[[nodiscard]] VelocityReconstruction prepare_velocity(const CubedSphereGrid& grid,
                                                      const std::vector<Vec3>& velocity,
                                                      const LimiterKind limiter) {
  VelocityReconstruction result{
      .gradient = least_squares_vector_gradient(grid, velocity),
      .factor = std::vector<Real>(velocity.size(), 1.0)};
  if (limiter == LimiterKind::kNone) return result;
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
    Real maximum_speed = norm(velocity[cell]);
    for (const auto& edge : grid.cell_cache()[cell].edges)
      maximum_speed = std::max(maximum_speed, norm(velocity[edge.neighbor]));
    for (const auto& cached_edge : grid.cell_cache()[cell].edges) {
      const auto& edge = grid.edges()[cached_edge.edge];
      const auto& edge_cache = grid.edge_cache()[cached_edge.edge];
      const EdgeTangentBasis basis{edge_cache.normal, edge_cache.tangent};
      const Vec3 center_at_face = project_tangent(velocity[cell], edge.center);
      const Vec3 increment =
          reconstruct_tangent_vector_cached(grid, cell, velocity[cell], edge.center,
                                            cached_edge.normalized_face_displacement_m,
                                            result.gradient[cell]) -
          center_at_face;
      Real minimum_normal = dot(center_at_face, basis.normal);
      Real maximum_normal = minimum_normal;
      Real minimum_tangent = dot(center_at_face, basis.tangent);
      Real maximum_tangent = minimum_tangent;
      for (const auto& neighbor_edge : grid.cell_cache()[cell].edges) {
        const Vec3 at_face =
            project_tangent(velocity[neighbor_edge.neighbor], edge.center);
        minimum_normal = std::min(minimum_normal, dot(at_face, basis.normal));
        maximum_normal = std::max(maximum_normal, dot(at_face, basis.normal));
        minimum_tangent = std::min(minimum_tangent, dot(at_face, basis.tangent));
        maximum_tangent = std::max(maximum_tangent, dot(at_face, basis.tangent));
      }
      result.factor[cell] =
          std::min(result.factor[cell], barth_factor(dot(center_at_face, basis.normal),
                                                     dot(increment, basis.normal),
                                                     minimum_normal, maximum_normal));
      result.factor[cell] =
          std::min(result.factor[cell], barth_factor(dot(center_at_face, basis.tangent),
                                                     dot(increment, basis.tangent),
                                                     minimum_tangent, maximum_tangent));
      result.factor[cell] = std::min(
          result.factor[cell], speed_factor(center_at_face, increment, maximum_speed));
    }
    result.factor[cell] = std::clamp(result.factor[cell], 0.0, 1.0);
  }
  return result;
}

[[nodiscard]] Vec3 reconstruct_velocity(const CubedSphereGrid& grid,
                                        const std::size_t cell,
                                        const EdgeGeometry& edge,
                                        const CachedCellEdgeGeometry& cached_edge,
                                        const std::vector<Vec3>& velocity,
                                        const VelocityReconstruction& prepared) {
  const Vec3 center_at_face = project_tangent(velocity[cell], edge.center);
  const Vec3 unlimited = reconstruct_tangent_vector_cached(
      grid, cell, velocity[cell], edge.center,
      cached_edge.normalized_face_displacement_m, prepared.gradient[cell]);
  return center_at_face + prepared.factor[cell] * (unlimited - center_at_face);
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
        const auto& cached = grid.edge_cache()[edge.id];
        const auto left = cached.left_cell;
        const auto right = cached.right_cell;
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
    const auto velocity_reconstruction = prepare_velocity(grid, velocity, limiter);
    for (std::size_t cell = 0; cell < cells; ++cell) {
      if (mass_reconstruction.factor[cell] < 1.0 - 1.0e-14 ||
          theta_reconstruction.factor[cell] < 1.0 - 1.0e-14 ||
          tracer_reconstruction.factor[cell] < 1.0 - 1.0e-14 ||
          temperature_reconstruction.factor[cell] < 1.0 - 1.0e-14 ||
          velocity_reconstruction.factor[cell] < 1.0 - 1.0e-14)
        ++result.limiter_activations;
    }
    for (const auto& edge : grid.edges()) {
      const auto& cached_edge = grid.edge_cache()[edge.id];
      const auto left = cached_edge.left_cell;
      const auto right = cached_edge.right_cell;
      const auto face = [&](const std::size_t cell) {
        const auto& cell_edge = cached_cell_edge(grid, cell, edge.id);
        return DryHydrostaticPrimitive{
            .air_mass_kg_m2 =
                reconstruct_scalar(cell, cell_edge, mass, mass_reconstruction),
            .velocity_m_s = reconstruct_velocity(grid, cell, edge, cell_edge, velocity,
                                                 velocity_reconstruction),
            .potential_temperature_k =
                reconstruct_scalar(cell, cell_edge, theta, theta_reconstruction),
            .tracer_mixing_ratio =
                reconstruct_scalar(cell, cell_edge, tracer, tracer_reconstruction),
            .temperature_k = reconstruct_scalar(cell, cell_edge, temperature,
                                                temperature_reconstruction)};
      };
      result.edge_levels[edge.id * levels + level] = {.left = face(left),
                                                      .right = face(right)};
    }
  }
  return result;
}

}  // namespace mps
