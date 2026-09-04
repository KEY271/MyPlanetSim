#include "myplanetsim/diagnostics/dry_hydrostatic_diagnostics.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <ostream>
#include <stdexcept>

#include "myplanetsim/numerics/spherical_operators.hpp"
namespace mps {
DryHydrostaticDiagnostics diagnose_dry_hydrostatic_budgets(
    const CubedSphereGrid& grid, const DryHydrostaticState& s,
    const DryHydrostaticDerived& d, const PlanetParameters& p) {
  DryHydrostaticDiagnostics x{0,
                              0,
                              0,
                              0,
                              0,
                              std::numeric_limits<Real>::infinity(),
                              -std::numeric_limits<Real>::infinity(),
                              std::numeric_limits<Real>::infinity(),
                              -std::numeric_limits<Real>::infinity()};
  const Real cv = p.heat_capacity_cp_j_kg_k - p.gas_constant_j_kg_k;
  for (std::size_t c = 0; c < d.cells; ++c) {
    auto area = grid.cells()[c].area_m2;
    auto pos = p.radius_m * grid.cells()[c].center;
    x.minimum_surface_pressure_pa =
        std::min(x.minimum_surface_pressure_pa, s.surface_pressure_pa[c]);
    x.maximum_surface_pressure_pa =
        std::max(x.maximum_surface_pressure_pa, s.surface_pressure_pa[c]);
    for (std::size_t k = 0; k < d.levels; ++k) {
      auto n = dry_hydrostatic_offset(c, k, d.levels);
      auto mass = area * d.air_mass_kg_m2[n];
      x.dry_mass_kg += mass;
      x.potential_temperature_mass_k_kg +=
          area * s.potential_temperature_mass_k_kg_m2[n];
      x.tracer_mass_kg += area * s.tracer_mass_kg_m2[n];
      x.total_energy_j += mass * (.5 * norm_squared(d.velocity_m_s[n]) +
                                  cv * d.temperature_k[n] + d.geopotential_m2_s2[n]);
      x.axial_angular_momentum_kg_m2_s +=
          mass * (cross(pos, d.velocity_m_s[n]).z +
                  p.rotation_rate_rad_s * (pos.x * pos.x + pos.y * pos.y));
      x.minimum_temperature_k = std::min(x.minimum_temperature_k, d.temperature_k[n]);
      x.maximum_temperature_k = std::max(x.maximum_temperature_k, d.temperature_k[n]);
    }
  }
  return x;
}

TerrainDiagnostics diagnose_terrain_budgets(
    const CubedSphereGrid& grid, const DryHydrostaticState& state,
    const DryHydrostaticDerived& d,
    const DryHydrostaticSources& sources,
    const std::span<const Real> surface_geopotential_m2_s2,
    const PlanetParameters& p) {
  const auto volume = d.cells * d.levels;
  if (d.cells != grid.cell_count() || d.levels == 0 ||
      state.surface_pressure_pa.size() != d.cells ||
      surface_geopotential_m2_s2.size() != d.cells ||
      sources.pressure_gradient_kg_m_s2.size() != volume)
    throw std::invalid_argument("terrain diagnostic shape mismatch");
  const auto surface_gradient =
      least_squares_gradient(grid, surface_geopotential_m2_s2);
  TerrainDiagnostics out{std::numeric_limits<Real>::infinity(),
                         -std::numeric_limits<Real>::infinity(),
                         0, 0, 0, 0, 0, 0, 0, 0, 0};
  Real weighted_l1 = 0;
  Real weighted_l2 = 0;
  for (std::size_t cell = 0; cell < d.cells; ++cell) {
    const Real area = grid.cells()[cell].area_m2;
    const Vec3 position = p.radius_m * grid.cells()[cell].center;
    const Real height = surface_geopotential_m2_s2[cell] / p.gravity_m_s2;
    out.minimum_surface_height_m = std::min(out.minimum_surface_height_m, height);
    out.maximum_surface_height_m = std::max(out.maximum_surface_height_m, height);
    out.maximum_surface_slope =
        std::max(out.maximum_surface_slope,
                 norm(surface_gradient[cell]) / p.gravity_m_s2);
    Vec3 column_force{};
    for (std::size_t level = 0; level < d.levels; ++level) {
      const auto offset = dry_hydrostatic_offset(cell, level, d.levels);
      const Real acceleration =
          norm(sources.pressure_gradient_kg_m_s2[offset]) /
          d.air_mass_kg_m2[offset];
      weighted_l1 += area * acceleration;
      weighted_l2 += area * acceleration * acceleration;
      out.pressure_gradient_linf_m_s2 =
          std::max(out.pressure_gradient_linf_m_s2, acceleration);
      out.maximum_wind_m_s =
          std::max(out.maximum_wind_m_s, norm(d.velocity_m_s[offset]));
      column_force = column_force + sources.pressure_gradient_kg_m_s2[offset];
      out.terrain_pressure_work_w +=
          area * dot(d.velocity_m_s[offset],
                     -d.air_mass_kg_m2[offset] * surface_gradient[cell]);
    }
    out.model_pressure_gradient_axial_torque_n_m +=
        area * cross(position, column_force).z;
    const Vec3 boundary_force =
        -(state.surface_pressure_pa[cell] / p.gravity_m_s2) *
        surface_gradient[cell];
    out.boundary_axial_torque_n_m +=
        area * cross(position, boundary_force).z;
  }
  const Real normalization = grid.total_area_m2() * static_cast<Real>(d.levels);
  out.pressure_gradient_l1_m_s2 = weighted_l1 / normalization;
  out.pressure_gradient_l2_m_s2 = std::sqrt(weighted_l2 / normalization);
  out.axial_torque_residual_n_m =
      out.model_pressure_gradient_axial_torque_n_m -
      out.boundary_axial_torque_n_m;
  return out;
}

Real absolute_pressure_velocity_pa_s(
    const Real b_half, const Real surface_pressure_tendency_pa_s,
    const Vec3 interface_velocity_m_s,
    const Vec3 interface_pressure_gradient_pa_m,
    const Real relative_mass_flux_kg_m2_s, const Real gravity_m_s2) {
  if (!std::isfinite(b_half) || !std::isfinite(surface_pressure_tendency_pa_s) ||
      !is_finite(interface_velocity_m_s) ||
      !is_finite(interface_pressure_gradient_pa_m) ||
      !std::isfinite(relative_mass_flux_kg_m2_s) || !(gravity_m_s2 > 0.0))
    throw std::invalid_argument("absolute pressure velocity arguments are invalid");
  return b_half * surface_pressure_tendency_pa_s +
         dot(interface_velocity_m_s, interface_pressure_gradient_pa_m) +
         gravity_m_s2 * relative_mass_flux_kg_m2_s;
}

void write_terrain_diagnostics(std::ostream& output,
                               const TerrainDiagnostics& diagnostics) {
  output << std::setprecision(std::numeric_limits<Real>::max_digits10)
         << "terrain.minimum_surface_height_m = "
         << diagnostics.minimum_surface_height_m << '\n'
         << "terrain.maximum_surface_height_m = "
         << diagnostics.maximum_surface_height_m << '\n'
         << "terrain.maximum_surface_slope = " << diagnostics.maximum_surface_slope
         << '\n'
         << "terrain.pressure_gradient_l1_m_s2 = "
         << diagnostics.pressure_gradient_l1_m_s2 << '\n'
         << "terrain.pressure_gradient_l2_m_s2 = "
         << diagnostics.pressure_gradient_l2_m_s2 << '\n'
         << "terrain.pressure_gradient_linf_m_s2 = "
         << diagnostics.pressure_gradient_linf_m_s2 << '\n'
         << "terrain.maximum_wind_m_s = " << diagnostics.maximum_wind_m_s << '\n'
         << "terrain.model_pressure_gradient_axial_torque_n_m = "
         << diagnostics.model_pressure_gradient_axial_torque_n_m << '\n'
         << "terrain.boundary_axial_torque_n_m = "
         << diagnostics.boundary_axial_torque_n_m << '\n'
         << "terrain.axial_torque_residual_n_m = "
         << diagnostics.axial_torque_residual_n_m << '\n'
         << "terrain.pressure_work_w = " << diagnostics.terrain_pressure_work_w
         << '\n';
  if (!output) throw std::runtime_error("failed while writing terrain diagnostics");
}
}  // namespace mps
