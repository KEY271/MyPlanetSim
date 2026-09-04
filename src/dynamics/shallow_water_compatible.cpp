#include "myplanetsim/dynamics/shallow_water_compatible.hpp"

#include <algorithm>
#include <cmath>
#include <span>
#include <stdexcept>
#include <vector>

#include "myplanetsim/core/validation.hpp"
#include "myplanetsim/dynamics/shallow_water_diffusion.hpp"
#include "myplanetsim/numerics/compatible_operators.hpp"

namespace mps {
namespace {

[[nodiscard]] ShallowWaterTendency zero_tendency(const std::size_t cell_count) {
  return {.depth = std::vector<Real>(cell_count),
          .momentum = std::vector<Vec3>(cell_count)};
}

[[nodiscard]] Real interpolate_to_edge(const CubedSphereGrid& grid,
                                       const DualEdgeGeometry& metric,
                                       const EdgeGeometry& edge,
                                       const std::span<const Real> cell_values) {
  const Real span = metric.left_center_to_edge_m + metric.right_center_to_edge_m;
  return (metric.right_center_to_edge_m * cell_values[grid.cell_index(edge.left_cell)] +
          metric.left_center_to_edge_m *
              cell_values[grid.cell_index(edge.right_cell)]) /
         span;
}

// Edge-normal derivative of a cell field. The difference along the dual edge
// only gives the derivative along the skewed dual direction, so the tangential
// derivative taken between the two primal vertices removes the skewness.
[[nodiscard]] std::vector<Real> edge_normal_derivative(
    const CubedSphereGrid& grid, const CubedSphereDualTopology& dual,
    const std::span<const CompatibleEdgeMetric> metrics,
    const std::span<const Real> cell_values) {
  const auto vertex_values = interpolate_cells_to_vertices(grid, dual, cell_values);
  std::vector<Real> derivative(grid.edge_count());
  for (const auto& edge : grid.edges()) {
    const Real along_dual = (cell_values[grid.cell_index(edge.right_cell)] -
                             cell_values[grid.cell_index(edge.left_cell)]) /
                            dual.edges()[edge.id].dual_length_m;
    // The primal edge runs from its first to its second vertex, which is the
    // edge tangent basis direction.
    const Real along_tangent =
        (vertex_values[edge.second_vertex] - vertex_values[edge.first_vertex]) /
        edge.length_m;
    derivative[edge.id] =
        (along_dual - metrics[edge.id].tangent_projection * along_tangent) /
        metrics[edge.id].normal_projection;
  }
  return derivative;
}

// Cartesian momentum tendency of a normal-velocity tendency, including the term
// that keeps the momentum consistent with the depth tendency.
void accumulate_momentum(const CubedSphereGrid& grid, const ShallowWaterState& state,
                         const std::span<const Real> edge_acceleration,
                         ShallowWaterTendency& tendency) {
  const auto acceleration = perot_cell_velocity(grid, edge_acceleration);
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
    tendency.momentum[cell] = project_tangent(
        tendency.momentum[cell] + state.depth[cell] * acceleration[cell],
        grid.cells()[cell].center);
  }
}

}  // namespace

ShallowWaterRhsComponents assemble_compatible_shallow_water_rhs(
    const CubedSphereGrid& grid, const CubedSphereDualTopology& dual,
    const ShallowWaterState& state, const ShallowWaterParameters& parameters,
    const Real gravity_m_s2, const Vec3 rotation_vector_rad_s) {
  require_positive(gravity_m_s2, "shallow-water gravity");
  if (!is_finite(rotation_vector_rad_s)) {
    throw std::invalid_argument("shallow-water rotation vector is non-finite");
  }
  validate_shallow_water_state(grid, state, parameters.depth_floor_m);
  std::vector<Vec3> velocity(grid.cell_count());
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
    velocity[cell] = state.velocity(cell);
  }
  const auto normal_velocity = edge_normal_velocity(grid, dual, velocity);
  const auto kinetic_energy = cell_kinetic_energy(grid, normal_velocity);
  const auto potential_vorticity = edge_potential_vorticity(
      grid, vertex_potential_vorticity(grid, dual, state.depth, normal_velocity,
                                       rotation_vector_rad_s));

  std::vector<Real> mass_flux(grid.edge_count());
  for (const auto& edge : grid.edges()) {
    const Real depth =
        interpolate_to_edge(grid, dual.edges()[edge.id], edge, state.depth);
    mass_flux[edge.id] = depth * normal_velocity[edge.id];
  }
  // The tangential mass flux of the potential vorticity flux term.
  const auto cell_mass_flux = perot_cell_velocity(grid, mass_flux);

  ShallowWaterRhsComponents result{
      .flux = zero_tendency(grid.cell_count()),
      .pressure = zero_tendency(grid.cell_count()),
      .coriolis = zero_tendency(grid.cell_count()),
      .diffusion = zero_tendency(grid.cell_count()),
      .orography = zero_tendency(grid.cell_count()),
      .total = zero_tendency(grid.cell_count()),
      .maximum_wave_speed_m_s = 0.0,
      .limiter_activations = 0,
  };
  const auto metrics = compatible_edge_metrics(grid, dual);
  std::vector<Real> free_surface(grid.cell_count());
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
    free_surface[cell] = gravity_m_s2 * state.depth[cell];
  }
  auto kinetic_gradient = edge_normal_derivative(grid, dual, metrics, kinetic_energy);
  auto pressure_gradient = edge_normal_derivative(grid, dual, metrics, free_surface);
  for (std::size_t edge = 0; edge < grid.edge_count(); ++edge) {
    kinetic_gradient[edge] = -kinetic_gradient[edge];
    pressure_gradient[edge] = -pressure_gradient[edge];
  }
  std::vector<Real> vorticity_flux(grid.edge_count());
  for (const auto& edge : grid.edges()) {
    const auto& metric = dual.edges()[edge.id];
    const std::size_t left = grid.cell_index(edge.left_cell);
    const std::size_t right = grid.cell_index(edge.right_cell);
    const Real integrated_mass = edge.length_m * mass_flux[edge.id];
    result.flux.depth[left] -= integrated_mass;
    result.flux.depth[right] += integrated_mass;

    const auto basis = edge_tangent_basis(edge);
    const Real span = metric.left_center_to_edge_m + metric.right_center_to_edge_m;
    const Vec3 flux_vector = (metric.right_center_to_edge_m * cell_mass_flux[left] +
                              metric.left_center_to_edge_m * cell_mass_flux[right]) /
                             span;
    vorticity_flux[edge.id] =
        potential_vorticity[edge.id] * dot(flux_vector, basis.tangent);
    result.maximum_wave_speed_m_s = std::max(
        result.maximum_wave_speed_m_s,
        std::abs(normal_velocity[edge.id]) +
            std::sqrt(gravity_m_s2 * std::max(state.depth[left], state.depth[right])));
  }
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
    result.flux.depth[cell] /= grid.cells()[cell].area_m2;
    result.flux.momentum[cell] = result.flux.depth[cell] * velocity[cell];
  }
  accumulate_momentum(grid, state, kinetic_gradient, result.flux);
  accumulate_momentum(grid, state, pressure_gradient, result.pressure);
  accumulate_momentum(grid, state, vorticity_flux, result.coriolis);
  result.diffusion = shallow_water_diffusion_tendency(
      grid, state, parameters.diffusion_kind, parameters.diffusion_coefficient);
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
    result.total.depth[cell] = result.flux.depth[cell] + result.diffusion.depth[cell];
    result.total.momentum[cell] = project_tangent(
        result.flux.momentum[cell] + result.pressure.momentum[cell] +
            result.coriolis.momentum[cell] + result.diffusion.momentum[cell],
        grid.cells()[cell].center);
  }
  return result;
}

}  // namespace mps
