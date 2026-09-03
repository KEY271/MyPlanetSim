#include "myplanetsim/diagnostics/shallow_water_diagnostics.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "myplanetsim/core/validation.hpp"
#include "myplanetsim/diagnostics/reductions.hpp"
#include "myplanetsim/numerics/spherical_operators.hpp"

namespace mps::diagnostics {
namespace {

void validate_tendency(const CubedSphereGrid& grid,
                       const ShallowWaterTendency& tendency) {
  if (tendency.depth.size() != grid.cell_count() ||
      tendency.momentum.size() != grid.cell_count()) {
    throw std::invalid_argument("shallow-water tendency size does not match grid");
  }
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
    require_finite(tendency.depth[cell], "depth tendency");
    if (!is_finite(tendency.momentum[cell])) {
      throw std::invalid_argument("momentum tendency is non-finite");
    }
  }
}

[[nodiscard]] ShallowWaterInvariantRates subtract_rates(
    const ShallowWaterInvariantRates first, const ShallowWaterInvariantRates second) {
  return {.mass = first.mass - second.mass,
          .energy = first.energy - second.energy,
          .potential_enstrophy = first.potential_enstrophy - second.potential_enstrophy,
          .axial_angular_momentum =
              first.axial_angular_momentum - second.axial_angular_momentum};
}

[[nodiscard]] ShallowWaterInvariantRates add_rates(
    const ShallowWaterInvariantRates first, const ShallowWaterInvariantRates second) {
  return {.mass = first.mass + second.mass,
          .energy = first.energy + second.energy,
          .potential_enstrophy = first.potential_enstrophy + second.potential_enstrophy,
          .axial_angular_momentum =
              first.axial_angular_momentum + second.axial_angular_momentum};
}

}  // namespace

ShallowWaterInvariants diagnose_shallow_water(const CubedSphereGrid& grid,
                                              const ShallowWaterState& state,
                                              const Real gravity_m_s2,
                                              const Vec3 rotation_vector_rad_s) {
  require_positive(gravity_m_s2, "shallow-water gravity");
  if (!is_finite(rotation_vector_rad_s)) {
    throw std::invalid_argument("shallow-water rotation vector is non-finite");
  }
  validate_shallow_water_state(grid, state, 0.0);
  std::vector<Vec3> velocity(grid.cell_count());
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
    velocity[cell] = state.velocity(cell);
  }
  const auto pv = shallow_water_potential_vorticity(grid, state.depth, velocity,
                                                    rotation_vector_rad_s);
  const Real rotation_rate_rad_s = norm(rotation_vector_rad_s);
  const Vec3 rotation_axis = rotation_rate_rad_s > 0.0
                                 ? rotation_vector_rad_s / rotation_rate_rad_s
                                 : Vec3{0.0, 0.0, 1.0};
  std::vector<Real> mass(grid.cell_count());
  std::vector<Real> energy(grid.cell_count());
  std::vector<Real> enstrophy(grid.cell_count());
  std::vector<Real> angular_momentum(grid.cell_count());
  std::vector<Real> radial_squared(grid.cell_count());
  Real maximum_radial = 0.0;
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
    const auto& geometry = grid.cells()[cell];
    const Real area = geometry.area_m2;
    const Real depth = state.depth[cell];
    const Real radial = dot(state.momentum[cell], geometry.center);
    mass[cell] = area * depth;
    energy[cell] = area * (0.5 * depth * norm_squared(velocity[cell]) +
                           0.5 * gravity_m_s2 * depth * depth);
    enstrophy[cell] = 0.5 * area * depth * pv[cell] * pv[cell];
    const Real axial_coordinate = dot(rotation_axis, geometry.center);
    const Real planetary = rotation_rate_rad_s * grid.radius_m() * grid.radius_m() *
                           (1.0 - axial_coordinate * axial_coordinate);
    angular_momentum[cell] =
        area * (grid.radius_m() *
                    dot(cross(geometry.center, state.momentum[cell]), rotation_axis) +
                depth * planetary);
    radial_squared[cell] = area * radial * radial;
    maximum_radial = std::max(maximum_radial, std::abs(radial));
  }
  const auto depth_range = min_max(state.depth);
  const auto pv_range = min_max(pv);
  return {.mass = compensated_sum(mass),
          .energy = compensated_sum(energy),
          .potential_enstrophy = compensated_sum(enstrophy),
          .axial_angular_momentum = compensated_sum(angular_momentum),
          .minimum_depth = depth_range.minimum,
          .maximum_depth = depth_range.maximum,
          .minimum_pv = pv_range.minimum,
          .maximum_pv = pv_range.maximum,
          .maximum_radial_momentum = maximum_radial,
          .rms_radial_momentum =
              std::sqrt(compensated_sum(radial_squared) / grid.total_area_m2())};
}

ShallowWaterInvariantRates shallow_water_invariant_rates(
    const CubedSphereGrid& grid, const ShallowWaterState& state,
    const ShallowWaterTendency& tendency, const Real gravity_m_s2,
    const Vec3 rotation_vector_rad_s) {
  require_positive(gravity_m_s2, "shallow-water gravity");
  if (!is_finite(rotation_vector_rad_s)) {
    throw std::invalid_argument("shallow-water rotation vector is non-finite");
  }
  validate_shallow_water_state(grid, state, 0.0);
  validate_tendency(grid, tendency);
  std::vector<Vec3> velocity(grid.cell_count());
  std::vector<Vec3> velocity_tendency(grid.cell_count());
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
    velocity[cell] = state.velocity(cell);
    velocity_tendency[cell] =
        (tendency.momentum[cell] - velocity[cell] * tendency.depth[cell]) /
        state.depth[cell];
  }
  const auto vorticity = finite_volume_curl(grid, velocity);
  const auto vorticity_tendency = finite_volume_curl(grid, velocity_tendency);
  const Real rotation_rate_rad_s = norm(rotation_vector_rad_s);
  const Vec3 rotation_axis = rotation_rate_rad_s > 0.0
                                 ? rotation_vector_rad_s / rotation_rate_rad_s
                                 : Vec3{0.0, 0.0, 1.0};
  std::vector<Real> mass(grid.cell_count());
  std::vector<Real> energy(grid.cell_count());
  std::vector<Real> enstrophy(grid.cell_count());
  std::vector<Real> angular_momentum(grid.cell_count());
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
    const auto& geometry = grid.cells()[cell];
    const Real area = geometry.area_m2;
    const Real depth = state.depth[cell];
    const Real depth_rate = tendency.depth[cell];
    const Real pv =
        (vorticity[cell] + 2.0 * dot(rotation_vector_rad_s, geometry.center)) / depth;
    mass[cell] = area * depth_rate;
    energy[cell] = area * ((gravity_m_s2 * depth - 0.5 * norm_squared(velocity[cell])) *
                               depth_rate +
                           dot(velocity[cell], tendency.momentum[cell]));
    enstrophy[cell] =
        area * (-0.5 * pv * pv * depth_rate + pv * vorticity_tendency[cell]);
    const Real axial_coordinate = dot(rotation_axis, geometry.center);
    const Real planetary = rotation_rate_rad_s * grid.radius_m() * grid.radius_m() *
                           (1.0 - axial_coordinate * axial_coordinate);
    angular_momentum[cell] =
        area * (grid.radius_m() * dot(cross(geometry.center, tendency.momentum[cell]),
                                      rotation_axis) +
                depth_rate * planetary);
  }
  return {.mass = compensated_sum(mass),
          .energy = compensated_sum(energy),
          .potential_enstrophy = compensated_sum(enstrophy),
          .axial_angular_momentum = compensated_sum(angular_momentum)};
}

ShallowWaterBudget make_shallow_water_budget(
    const CubedSphereGrid& grid, const ShallowWaterState& state,
    const ShallowWaterTendency& flux, const ShallowWaterTendency& coriolis,
    const ShallowWaterTendency& pressure, const ShallowWaterTendency& diffusion,
    const Real gravity_m_s2, const Vec3 rotation_vector_rad_s) {
  const auto flux_rates = shallow_water_invariant_rates(grid, state, flux, gravity_m_s2,
                                                        rotation_vector_rad_s);
  const auto coriolis_rates = shallow_water_invariant_rates(
      grid, state, coriolis, gravity_m_s2, rotation_vector_rad_s);
  const auto pressure_rates = shallow_water_invariant_rates(
      grid, state, pressure, gravity_m_s2, rotation_vector_rad_s);
  const auto diffusion_rates = shallow_water_invariant_rates(
      grid, state, diffusion, gravity_m_s2, rotation_vector_rad_s);
  ShallowWaterTendency total{.depth = std::vector<Real>(grid.cell_count()),
                             .momentum = std::vector<Vec3>(grid.cell_count())};
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
    total.depth[cell] = flux.depth[cell] + coriolis.depth[cell] + pressure.depth[cell] +
                        diffusion.depth[cell];
    total.momentum[cell] = flux.momentum[cell] + coriolis.momentum[cell] +
                           pressure.momentum[cell] + diffusion.momentum[cell];
  }
  const auto total_rates = shallow_water_invariant_rates(
      grid, state, total, gravity_m_s2, rotation_vector_rad_s);
  const auto component_sum = add_rates(add_rates(flux_rates, coriolis_rates),
                                       add_rates(pressure_rates, diffusion_rates));
  return {.flux = flux_rates,
          .coriolis = coriolis_rates,
          .pressure = pressure_rates,
          .diffusion = diffusion_rates,
          .total = total_rates,
          .residual = subtract_rates(total_rates, component_sum)};
}

}  // namespace mps::diagnostics
