#include "myplanetsim/dynamics/shallow_water_rhs.hpp"

#include <algorithm>
#include <stdexcept>
#include <vector>

#include "myplanetsim/core/validation.hpp"

namespace mps {
namespace {

[[nodiscard]] ShallowWaterTendency zero_tendency(const std::size_t cell_count) {
  return {.depth = std::vector<Real>(cell_count),
          .momentum = std::vector<Vec3>(cell_count)};
}

void scatter_edge(const CubedSphereGrid& grid, const EdgeGeometry& edge,
                  const Real mass_flux, const Vec3 momentum_flux,
                  ShallowWaterTendency& tendency) {
  const std::size_t left = grid.cell_index(edge.left_cell);
  const std::size_t right = grid.cell_index(edge.right_cell);
  const Real integrated_mass = edge.length_m * mass_flux;
  const Vec3 integrated_momentum = edge.length_m * momentum_flux;
  tendency.depth[left] -= integrated_mass;
  tendency.depth[right] += integrated_mass;
  tendency.momentum[left] = tendency.momentum[left] - integrated_momentum;
  tendency.momentum[right] = tendency.momentum[right] + integrated_momentum;
}

void divide_by_area_and_project(const CubedSphereGrid& grid,
                                ShallowWaterTendency& tendency) {
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
    const Real inverse_area = 1.0 / grid.cells()[cell].area_m2;
    tendency.depth[cell] *= inverse_area;
    tendency.momentum[cell] = project_tangent(inverse_area * tendency.momentum[cell],
                                              grid.cells()[cell].center);
  }
}

}  // namespace

ShallowWaterRhsComponents assemble_first_order_shallow_water_rhs(
    const CubedSphereGrid& grid, const ShallowWaterState& state,
    const Real gravity_m_s2, const Real rotation_rate_rad_s, const Real depth_floor_m) {
  require_positive(gravity_m_s2, "shallow-water gravity");
  require_finite(rotation_rate_rad_s, "shallow-water rotation rate");
  validate_shallow_water_state(grid, state, depth_floor_m);
  ShallowWaterRhsComponents result{
      .flux = zero_tendency(grid.cell_count()),
      .pressure = zero_tendency(grid.cell_count()),
      .coriolis = zero_tendency(grid.cell_count()),
      .diffusion = zero_tendency(grid.cell_count()),
      .total = zero_tendency(grid.cell_count()),
      .maximum_wave_speed_m_s = 0.0,
  };
  for (const auto& edge : grid.edges()) {
    const std::size_t left = grid.cell_index(edge.left_cell);
    const std::size_t right = grid.cell_index(edge.right_cell);
    const auto basis = edge_tangent_basis(edge);
    const ShallowWaterPrimitive left_state{
        .depth_m = state.depth[left],
        .velocity_m_s = project_tangent(state.velocity(left), edge.center),
    };
    const ShallowWaterPrimitive right_state{
        .depth_m = state.depth[right],
        .velocity_m_s = project_tangent(state.velocity(right), edge.center),
    };
    const auto edge_flux =
        rusanov_shallow_water_flux(left_state, right_state, basis, gravity_m_s2);
    const Vec3 pressure_flux = 0.25 * gravity_m_s2 *
                               (left_state.depth_m * left_state.depth_m +
                                right_state.depth_m * right_state.depth_m) *
                               basis.normal;
    scatter_edge(grid, edge, edge_flux.mass_m2_s,
                 edge_flux.momentum_m3_s2 - pressure_flux, result.flux);
    scatter_edge(grid, edge, 0.0, pressure_flux, result.pressure);
    result.maximum_wave_speed_m_s =
        std::max(result.maximum_wave_speed_m_s, edge_flux.maximum_wave_speed_m_s);
  }
  divide_by_area_and_project(grid, result.flux);
  divide_by_area_and_project(grid, result.pressure);
  const Vec3 omega{0.0, 0.0, rotation_rate_rad_s};
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
    Vec3 pressure_geometry_correction{};
    for (const std::size_t edge_id : grid.cell_edges(grid.cell_id(cell))) {
      const auto& edge = grid.edge(edge_id);
      pressure_geometry_correction =
          pressure_geometry_correction +
          static_cast<Real>(grid.edge_sign_for_cell(edge_id, grid.cell_id(cell))) *
              edge.length_m * edge_tangent_basis(edge).normal;
    }
    const Real cell_pressure =
        0.5 * gravity_m_s2 * state.depth[cell] * state.depth[cell];
    result.pressure.momentum[cell] = project_tangent(
        result.pressure.momentum[cell] +
            (cell_pressure / grid.cells()[cell].area_m2) * pressure_geometry_correction,
        grid.cells()[cell].center);
    result.coriolis.momentum[cell] =
        -2.0 *
        project_tangent(cross(omega, state.momentum[cell]), grid.cells()[cell].center);
    result.total.depth[cell] = result.flux.depth[cell];
    result.total.momentum[cell] =
        project_tangent(result.flux.momentum[cell] + result.pressure.momentum[cell] +
                            result.coriolis.momentum[cell],
                        grid.cells()[cell].center);
  }
  return result;
}

}  // namespace mps
