#include "myplanetsim/dynamics/dry_hydrostatic_reconstruction.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#if defined(MPS_ENABLE_OPENMP)
#include <omp.h>
#endif

#include "myplanetsim/numerics/spherical_operators.hpp"

namespace mps {
namespace {

#if defined(MPS_ENABLE_OPENMP)
constexpr int kMaximumReconstructionThreads = 4;
#endif

[[nodiscard]] int reconstruction_thread_count(const std::size_t levels) noexcept {
#if defined(MPS_ENABLE_OPENMP)
  const int available = std::min(kMaximumReconstructionThreads, omp_get_max_threads());
  return static_cast<int>(
      std::min(levels, static_cast<std::size_t>(std::max(1, available))));
#else
  static_cast<void>(levels);
  return 1;
#endif
}

[[nodiscard]] int reconstruction_thread_index() noexcept {
#if defined(MPS_ENABLE_OPENMP)
  return omp_get_thread_num();
#else
  return 0;
#endif
}

struct ScalarReconstruction {
  std::span<Vec3> gradient;
  std::span<Real> factor;
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
  const Real maximum_squared = maximum_speed * maximum_speed;
  if (norm_squared(center + increment) <= maximum_squared) return 1.0;
  const Real a = norm_squared(increment);
  const Real b = 2.0 * dot(center, increment);
  const Real c = norm_squared(center) - maximum_squared;
  const Real discriminant = std::max(0.0, b * b - 4.0 * a * c);
  return std::clamp((-b + std::sqrt(discriminant)) / (2.0 * a), 0.0, 1.0);
}

[[nodiscard]] ScalarReconstruction prepare_scalar(const CubedSphereGrid& grid,
                                                  const std::span<const Real> values,
                                                  const LimiterKind limiter,
                                                  const std::span<Vec3> gradient,
                                                  const std::span<Real> factor) {
  least_squares_gradient(grid, values, gradient);
  std::fill(factor.begin(), factor.end(), 1.0);
  ScalarReconstruction result{.gradient = gradient, .factor = factor};
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
                                      const std::span<const Real> values,
                                      const ScalarReconstruction& prepared) {
  return values[cell] +
         prepared.factor[cell] * dot(prepared.gradient[cell], edge.face_displacement_m);
}

struct VelocityReconstruction {
  std::span<TangentVectorGradient> gradient;
  std::span<Real> factor;
};

[[nodiscard]] Vec3 reconstruct_velocity_unchecked(
    const CubedSphereGrid& grid, const std::size_t cell, const Vec3 cell_value,
    const Vec3 face, const Vec3 displacement, const TangentVectorGradient& gradient) {
  const auto& cached = grid.cell_cache()[cell];
  const Vec3 increment =
      gradient.alpha_derivative * dot(displacement, cached.basis.alpha) +
      gradient.beta_derivative * dot(displacement, cached.basis.beta);
  const Vec3 at_cell = cell_value + increment;
  const Vec3 center = grid.cells()[cell].center;
  return at_cell - (dot(at_cell, face) / (1.0 + dot(center, face))) * (center + face);
}

// Limits the reconstructed tangent velocity the same way the shallow-water path does:
// the face normal and tangent components stay inside the neighbourhood range and the
// face speed stays inside the neighbourhood maximum. Without this the unlimited
// least-squares vector gradient overshoots at cubed-sphere seams.
[[nodiscard]] VelocityReconstruction prepare_velocity(
    const CubedSphereGrid& grid, const std::span<const Vec3> velocity,
    const LimiterKind limiter, const std::span<TangentVectorGradient> gradient,
    const std::span<Real> factor) {
  least_squares_vector_gradient(grid, velocity, gradient);
  std::fill(factor.begin(), factor.end(), 1.0);
  VelocityReconstruction result{.gradient = gradient, .factor = factor};
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
          reconstruct_velocity_unchecked(grid, cell, velocity[cell], edge.center,
                                         cached_edge.normalized_face_displacement_m,
                                         result.gradient[cell]) -
          center_at_face;
      Real minimum_normal = dot(center_at_face, basis.normal);
      Real maximum_normal = minimum_normal;
      Real minimum_tangent = dot(center_at_face, basis.tangent);
      Real maximum_tangent = minimum_tangent;
      for (const auto& neighbor_edge : grid.cell_cache()[cell].edges) {
        const Vec3 neighbor_velocity = velocity[neighbor_edge.neighbor];
        minimum_normal = std::min(minimum_normal, dot(neighbor_velocity, basis.normal));
        maximum_normal = std::max(maximum_normal, dot(neighbor_velocity, basis.normal));
        minimum_tangent =
            std::min(minimum_tangent, dot(neighbor_velocity, basis.tangent));
        maximum_tangent =
            std::max(maximum_tangent, dot(neighbor_velocity, basis.tangent));
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
                                        const std::span<const Vec3> velocity,
                                        const VelocityReconstruction& prepared) {
  const Vec3 center_at_face = project_tangent(velocity[cell], edge.center);
  const Vec3 unlimited = reconstruct_velocity_unchecked(
      grid, cell, velocity[cell], edge.center,
      cached_edge.normalized_face_displacement_m, prepared.gradient[cell]);
  return center_at_face + prepared.factor[cell] * (unlimited - center_at_face);
}

void copy_level_scalar(const std::span<const Real> volume, const std::size_t cells,
                       const std::size_t levels, const std::size_t level,
                       const std::span<Real> result) {
  for (std::size_t cell = 0; cell < cells; ++cell)
    result[cell] = volume[dry_hydrostatic_offset(cell, level, levels)];
}

void copy_level_vector(const std::span<const Vec3> volume, const std::size_t cells,
                       const std::size_t levels, const std::size_t level,
                       const std::span<Vec3> result) {
  for (std::size_t cell = 0; cell < cells; ++cell)
    result[cell] = volume[dry_hydrostatic_offset(cell, level, levels)];
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
  DryHydrostaticReconstruction result;
  DryHydrostaticReconstructionWorkspace workspace;
  reconstruct_dry_hydrostatic_face_states(grid, derived, reconstruction, limiter,
                                          result, workspace);
  return result;
}

void reconstruct_dry_hydrostatic_face_states(
    const CubedSphereGrid& grid, const DryHydrostaticDerived& derived,
    const ReconstructionKind reconstruction, const LimiterKind limiter,
    DryHydrostaticReconstruction& result,
    DryHydrostaticReconstructionWorkspace& workspace) {
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

  result.levels = levels;
  result.edge_levels.resize(grid.edge_count() * levels);
  std::uint64_t limiter_activations = 0;
  const int thread_count = reconstruction_thread_count(levels);
  if (workspace.workers.size() < static_cast<std::size_t>(thread_count))
    workspace.workers.resize(static_cast<std::size_t>(thread_count));
  for (int thread = 0; thread < thread_count; ++thread) {
    auto& worker = workspace.workers[static_cast<std::size_t>(thread)];
    worker.mass.resize(cells);
    worker.velocity.resize(cells);
    worker.potential_temperature.resize(cells);
    worker.tracer.resize(cells);
    worker.temperature.resize(cells);
    for (auto& gradient : worker.scalar_gradients) gradient.resize(cells);
    for (auto& factor : worker.limiter_factors) factor.resize(cells);
    worker.velocity_gradient.resize(cells);
  }
#if defined(MPS_ENABLE_OPENMP)
#pragma omp parallel for schedule(static) num_threads(thread_count) \
    reduction(+ : limiter_activations)
#endif
  for (std::size_t level = 0; level < levels; ++level) {
    auto& worker =
        workspace.workers[static_cast<std::size_t>(reconstruction_thread_index())];
    copy_level_scalar(derived.air_mass_kg_m2, cells, levels, level, worker.mass);
    copy_level_vector(derived.velocity_m_s, cells, levels, level, worker.velocity);
    copy_level_scalar(derived.potential_temperature_k, cells, levels, level,
                      worker.potential_temperature);
    copy_level_scalar(derived.tracer_mixing_ratio, cells, levels, level, worker.tracer);
    copy_level_scalar(derived.temperature_k, cells, levels, level, worker.temperature);
    const auto& mass = worker.mass;
    const auto& velocity = worker.velocity;
    const auto& theta = worker.potential_temperature;
    const auto& tracer = worker.tracer;
    const auto& temperature = worker.temperature;

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

    const auto mass_reconstruction = prepare_scalar(
        grid, mass, limiter, worker.scalar_gradients[0], worker.limiter_factors[0]);
    const auto theta_reconstruction = prepare_scalar(
        grid, theta, limiter, worker.scalar_gradients[1], worker.limiter_factors[1]);
    const auto tracer_reconstruction = prepare_scalar(
        grid, tracer, limiter, worker.scalar_gradients[2], worker.limiter_factors[2]);
    const auto temperature_reconstruction =
        prepare_scalar(grid, temperature, limiter, worker.scalar_gradients[3],
                       worker.limiter_factors[3]);
    const auto velocity_reconstruction = prepare_velocity(
        grid, velocity, limiter, worker.velocity_gradient, worker.limiter_factors[4]);
    for (std::size_t cell = 0; cell < cells; ++cell) {
      if (mass_reconstruction.factor[cell] < 1.0 - 1.0e-14 ||
          theta_reconstruction.factor[cell] < 1.0 - 1.0e-14 ||
          tracer_reconstruction.factor[cell] < 1.0 - 1.0e-14 ||
          temperature_reconstruction.factor[cell] < 1.0 - 1.0e-14 ||
          velocity_reconstruction.factor[cell] < 1.0 - 1.0e-14)
        ++limiter_activations;
    }
    for (const auto& edge : grid.edges()) {
      const auto& cached_edge = grid.edge_cache()[edge.id];
      const auto left = cached_edge.left_cell;
      const auto right = cached_edge.right_cell;
      const auto face = [&](const std::size_t cell) {
        const std::size_t slot =
            cell == left ? cached_edge.left_slot : cached_edge.right_slot;
        const auto& cell_edge = grid.cell_cache()[cell].edges[slot];
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
  result.limiter_activations = limiter_activations;
}

}  // namespace mps
