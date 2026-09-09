#include "myplanetsim/dynamics/dry_hydrostatic_fast_operator.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

#if defined(MPS_ENABLE_OPENMP)
#include <omp.h>
#endif

#include "myplanetsim/core/validation.hpp"
#include "myplanetsim/numerics/spherical_operators.hpp"

namespace mps {
namespace {

[[nodiscard]] std::size_t fast_worker_count() {
#if defined(MPS_ENABLE_OPENMP)
  return static_cast<std::size_t>(std::max(1, omp_get_max_threads()));
#else
  return 1;
#endif
}

[[nodiscard]] std::size_t fast_worker_index() {
#if defined(MPS_ENABLE_OPENMP)
  return static_cast<std::size_t>(omp_get_thread_num());
#else
  return 0;
#endif
}

struct ColumnDiagnostics {
  std::vector<Real> pressure_pa;
  std::vector<Real> geopotential_m2_s2;
};

[[nodiscard]] ColumnDiagnostics diagnose_column(
    const AtmosphericHybridCoordinate& coordinate, const PlanetParameters& planet,
    const Real surface_pressure_pa,
    const std::span<const Real> potential_temperature_mass) {
  const auto geometry = coordinate.geometry(
      surface_pressure_pa, planet.gravity_m_s2, planet.gas_constant_j_kg_k,
      planet.heat_capacity_cp_j_kg_k, planet.reference_pressure_pa);
  if (potential_temperature_mass.size() != coordinate.levels())
    throw std::invalid_argument("reference-linear column shape mismatch");
  std::vector<Real> theta(coordinate.levels());
  for (std::size_t level = 0; level < coordinate.levels(); ++level)
    theta[level] = potential_temperature_mass[level] / geometry.air_mass_kg_m2[level];
  const auto hydrostatic = integrate_hydrostatic_column(
      geometry, theta, planet.heat_capacity_cp_j_kg_k, planet.gravity_m_s2, 0.0);
  return {.pressure_pa = geometry.pressure_full_pa,
          .geopotential_m2_s2 = hydrostatic.geopotential_full_m2_s2};
}

void validate_fast_shape(const DryHydrostaticFastTendency& value,
                         const std::size_t cells, const std::size_t levels) {
  const auto volume = cells * levels;
  if (value.surface_pressure_pa_s.size() != cells ||
      value.tendency.air_mass.size() != volume ||
      value.tendency.momentum.size() != volume ||
      value.tendency.potential_temperature_mass.size() != volume ||
      value.tendency.tracer_mass.size() != volume)
    throw std::invalid_argument("dry fast tendency shape mismatch");
}

[[nodiscard]] DryHydrostaticFastTendency combine_fast_tendencies(
    const DryHydrostaticFastTendency& left, const DryHydrostaticFastTendency& right,
    const Real right_scale) {
  const auto cells = left.surface_pressure_pa_s.size();
  if (cells == 0 || right.surface_pressure_pa_s.size() != cells ||
      left.tendency.air_mass.empty() || left.tendency.air_mass.size() % cells != 0)
    throw std::invalid_argument("dry fast tendencies cannot be combined");
  const auto levels = left.tendency.air_mass.size() / cells;
  validate_fast_shape(left, cells, levels);
  validate_fast_shape(right, cells, levels);
  DryHydrostaticFastTendency result = left;
  for (std::size_t cell = 0; cell < cells; ++cell)
    result.surface_pressure_pa_s[cell] +=
        right_scale * right.surface_pressure_pa_s[cell];
  for (std::size_t n = 0; n < cells * levels; ++n) {
    result.tendency.air_mass[n] += right_scale * right.tendency.air_mass[n];
    result.tendency.momentum[n] =
        result.tendency.momentum[n] + right_scale * right.tendency.momentum[n];
    result.tendency.potential_temperature_mass[n] +=
        right_scale * right.tendency.potential_temperature_mass[n];
    result.tendency.tracer_mass[n] += right_scale * right.tendency.tracer_mass[n];
  }
  return result;
}

}  // namespace

DryHydrostaticFastOperator make_dry_hydrostatic_fast_operator(
    const AtmosphericHybridCoordinate& coordinate, const PlanetParameters& planet,
    const DryHydrostaticReferenceColumn& reference) {
  planet.validate();
  const auto levels = coordinate.levels();
  if (levels == 0 || reference.geometry.air_mass_kg_m2.size() != levels ||
      reference.potential_temperature_k.size() != levels ||
      reference.potential_temperature_mass_k_kg_m2.size() != levels)
    throw std::invalid_argument("dry reference column shape mismatch");

  DryHydrostaticFastOperator result{
      .levels = levels,
      .b_half = coordinate.coefficients().b_half,
      .reference_air_mass_kg_m2 = reference.geometry.air_mass_kg_m2,
      .reference_potential_temperature_k = reference.potential_temperature_k,
      .reference_interface_potential_temperature_k = std::vector<Real>(levels + 1),
      .reference_specific_volume_m3_kg = std::vector<Real>(levels),
      .pressure_from_surface_pressure = std::vector<Real>(levels),
      .geopotential_from_surface_pressure = std::vector<Real>(levels),
      .geopotential_from_potential_temperature_mass =
          std::vector<Real>(levels * levels)};

  result.reference_interface_potential_temperature_k.front() =
      reference.potential_temperature_k.front();
  result.reference_interface_potential_temperature_k.back() =
      reference.potential_temperature_k.back();
  for (std::size_t interface = 1; interface < levels; ++interface) {
    result.reference_interface_potential_temperature_k[interface] =
        0.5 * (reference.potential_temperature_k[interface - 1] +
               reference.potential_temperature_k[interface]);
  }
  for (std::size_t level = 0; level < levels; ++level) {
    result.reference_specific_volume_m3_kg[level] =
        planet.gas_constant_j_kg_k * reference.temperature_k /
        reference.geometry.pressure_full_pa[level];
  }

  const auto& reference_q = reference.potential_temperature_mass_k_kg_m2;
  const Real pressure_step = std::max(1.0e-3, 1.0e-4 * reference.surface_pressure_pa);
  ColumnDiagnostics lower;
  ColumnDiagnostics upper;
  Real denominator = 2.0 * pressure_step;
  try {
    lower = diagnose_column(coordinate, planet,
                            reference.surface_pressure_pa - pressure_step, reference_q);
    upper = diagnose_column(coordinate, planet,
                            reference.surface_pressure_pa + pressure_step, reference_q);
  } catch (const std::invalid_argument&) {
    const auto center =
        diagnose_column(coordinate, planet, reference.surface_pressure_pa, reference_q);
    try {
      upper =
          diagnose_column(coordinate, planet,
                          reference.surface_pressure_pa + pressure_step, reference_q);
      lower = center;
    } catch (const std::invalid_argument&) {
      lower =
          diagnose_column(coordinate, planet,
                          reference.surface_pressure_pa - pressure_step, reference_q);
      upper = center;
    }
    denominator = pressure_step;
  }
  for (std::size_t level = 0; level < levels; ++level) {
    result.pressure_from_surface_pressure[level] =
        (upper.pressure_pa[level] - lower.pressure_pa[level]) / denominator;
    result.geopotential_from_surface_pressure[level] =
        (upper.geopotential_m2_s2[level] - lower.geopotential_m2_s2[level]) /
        denominator;
  }

  for (std::size_t source = 0; source < levels; ++source) {
    auto lower_q = reference_q;
    auto upper_q = reference_q;
    const Real step = std::max(1.0, 1.0e-6 * std::abs(reference_q[source]));
    lower_q[source] -= step;
    upper_q[source] += step;
    lower = diagnose_column(coordinate, planet, reference.surface_pressure_pa, lower_q);
    upper = diagnose_column(coordinate, planet, reference.surface_pressure_pa, upper_q);
    for (std::size_t level = 0; level < levels; ++level) {
      result.geopotential_from_potential_temperature_mass[level * levels + source] =
          (upper.geopotential_m2_s2[level] - lower.geopotential_m2_s2[level]) /
          (2.0 * step);
    }
  }
  return result;
}

DryHydrostaticFastPerturbation make_dry_hydrostatic_fast_perturbation(
    const DryHydrostaticState& state, const DryHydrostaticReferenceColumn& reference) {
  const auto cells = state.surface_pressure_pa.size();
  const auto levels = reference.geometry.air_mass_kg_m2.size();
  const auto volume = cells * levels;
  if (cells == 0 || levels == 0 ||
      state.horizontal_momentum_mass_kg_m_s.size() != volume ||
      state.potential_temperature_mass_k_kg_m2.size() != volume)
    throw std::invalid_argument("dry state and fast reference shapes differ");
  DryHydrostaticFastPerturbation result{
      .surface_pressure_pa = std::vector<Real>(cells),
      .horizontal_momentum_mass_kg_m_s = state.horizontal_momentum_mass_kg_m_s,
      .potential_temperature_mass_k_kg_m2 = std::vector<Real>(volume)};
  for (std::size_t cell = 0; cell < cells; ++cell) {
    result.surface_pressure_pa[cell] =
        state.surface_pressure_pa[cell] - reference.surface_pressure_pa;
    for (std::size_t level = 0; level < levels; ++level) {
      const auto n = dry_hydrostatic_offset(cell, level, levels);
      result.potential_temperature_mass_k_kg_m2[n] =
          state.potential_temperature_mass_k_kg_m2[n] -
          reference.potential_temperature_mass_k_kg_m2[level];
    }
  }
  return result;
}

void apply_dry_hydrostatic_fast_operator(
    const CubedSphereGrid& grid, const PlanetParameters& planet,
    const DryHydrostaticFastOperator& op,
    const DryHydrostaticFastPerturbation& perturbation,
    DryHydrostaticFastTendency& result, DryHydrostaticFastOperatorWorkspace& workspace,
    const bool compute_scalar_from_momentum, const bool compute_momentum_from_scalar) {
  const auto cells = grid.cell_count();
  const auto levels = op.levels;
  const auto volume = cells * levels;
  if (cells == 0 || levels == 0 || perturbation.surface_pressure_pa.size() != cells ||
      perturbation.horizontal_momentum_mass_kg_m_s.size() != volume ||
      perturbation.potential_temperature_mass_k_kg_m2.size() != volume ||
      op.b_half.size() != levels + 1 || op.reference_air_mass_kg_m2.size() != levels ||
      op.reference_potential_temperature_k.size() != levels ||
      op.reference_interface_potential_temperature_k.size() != levels + 1 ||
      op.reference_specific_volume_m3_kg.size() != levels ||
      op.pressure_from_surface_pressure.size() != levels ||
      op.geopotential_from_surface_pressure.size() != levels ||
      op.geopotential_from_potential_temperature_mass.size() != levels * levels)
    throw std::invalid_argument("dry fast operator shape mismatch");

  if (compute_scalar_from_momentum) {
    result.surface_pressure_pa_s.assign(cells, 0.0);
    result.tendency.air_mass.assign(volume, 0.0);
    result.tendency.potential_temperature_mass.assign(volume, 0.0);
    workspace.horizontal_air_mass_tendency.assign(volume, 0.0);
    workspace.horizontal_potential_temperature_mass_tendency.assign(volume, 0.0);
    const auto edges = grid.edges();
    workspace.edge_integrated_mass_flux.resize(edges.size() * levels);
#if defined(MPS_ENABLE_OPENMP)
#pragma omp parallel for schedule(static)
#endif
    for (std::ptrdiff_t edge_index = 0;
         edge_index < static_cast<std::ptrdiff_t>(edges.size()); ++edge_index) {
      const auto edge_id = static_cast<std::size_t>(edge_index);
      const auto& edge = edges[edge_id];
      const auto& cached = grid.edge_cache()[edge_id];
      for (std::size_t level = 0; level < levels; ++level) {
        const auto left = dry_hydrostatic_offset(cached.left_cell, level, levels);
        const auto right = dry_hydrostatic_offset(cached.right_cell, level, levels);
        workspace.edge_integrated_mass_flux[edge_id * levels + level] =
            0.5 *
            dot(perturbation.horizontal_momentum_mass_kg_m_s[left] +
                    perturbation.horizontal_momentum_mass_kg_m_s[right],
                edge.outward_normal_from_left) *
            edge.length_m;
      }
    }
#if defined(MPS_ENABLE_OPENMP)
#pragma omp parallel for schedule(static)
#endif
    for (std::ptrdiff_t cell_index = 0; cell_index < static_cast<std::ptrdiff_t>(cells);
         ++cell_index) {
      const auto cell = static_cast<std::size_t>(cell_index);
      std::array<std::size_t, 4> cell_edges{};
      for (std::size_t side = 0; side < cell_edges.size(); ++side)
        cell_edges[side] = grid.cell_cache()[cell].edges[side].edge;
      std::ranges::sort(cell_edges);
      for (const auto edge : cell_edges) {
        const auto& cached = grid.edge_cache()[edge];
        const Real sign = cached.left_cell == cell ? -1.0 : 1.0;
        for (std::size_t level = 0; level < levels; ++level) {
          const auto n = dry_hydrostatic_offset(cell, level, levels);
          const Real mass_tendency =
              sign * workspace.edge_integrated_mass_flux[edge * levels + level] /
              grid.cells()[cell].area_m2;
          workspace.horizontal_air_mass_tendency[n] += mass_tendency;
          workspace.horizontal_potential_temperature_mass_tendency[n] +=
              op.reference_potential_temperature_k[level] * mass_tendency;
        }
      }
    }

    workspace.vertical_mass_flux_workers.resize(fast_worker_count());
    workspace.failures.resize(cells);
#if defined(MPS_ENABLE_OPENMP)
#pragma omp parallel for schedule(static)
#endif
    for (std::ptrdiff_t cell_index = 0; cell_index < static_cast<std::ptrdiff_t>(cells);
         ++cell_index) {
      const auto cell = static_cast<std::size_t>(cell_index);
      const auto begin = cell * levels;
      auto& vertical_mass_flux =
          workspace.vertical_mass_flux_workers[fast_worker_index()];
      workspace.failures[cell] = nullptr;
      try {
        diagnose_vertical_mass_flux(
            std::span<const Real>(workspace.horizontal_air_mass_tendency.data() + begin,
                                  levels),
            op.b_half, planet.gravity_m_s2, vertical_mass_flux);
        result.surface_pressure_pa_s[cell] =
            vertical_mass_flux.surface_pressure_tendency_pa_s;
        for (std::size_t level = 0; level < levels; ++level) {
          const auto n = begin + level;
          result.tendency.air_mass[n] =
              vertical_mass_flux.target_air_mass_tendency_kg_m2_s[level];
          result.tendency.potential_temperature_mass[n] =
              workspace.horizontal_potential_temperature_mass_tendency[n] +
              vertical_mass_flux.interface_flux_kg_m2_s[level] *
                  op.reference_interface_potential_temperature_k[level] -
              vertical_mass_flux.interface_flux_kg_m2_s[level + 1] *
                  op.reference_interface_potential_temperature_k[level + 1];
        }
      } catch (...) {
        workspace.failures[cell] = std::current_exception();
      }
    }
    for (std::size_t cell = 0; cell < cells; ++cell)
      if (workspace.failures[cell] != nullptr)
        std::rethrow_exception(workspace.failures[cell]);
  }

  if (compute_momentum_from_scalar) {
    result.tendency.momentum.assign(volume, {});
    const auto workers = std::min(levels, fast_worker_count());
    workspace.level_workers.resize(workers);
    for (auto& worker : workspace.level_workers) {
      worker.geopotential_perturbation.resize(cells);
      worker.geopotential_gradient.resize(cells);
    }
    workspace.failures.resize(levels);
#if defined(MPS_ENABLE_OPENMP)
#pragma omp parallel for schedule(static) num_threads(workers)
#endif
    for (std::ptrdiff_t level_index = 0;
         level_index < static_cast<std::ptrdiff_t>(levels); ++level_index) {
      const auto level = static_cast<std::size_t>(level_index);
      auto& worker = workspace.level_workers[fast_worker_index()];
      workspace.failures[level] = nullptr;
      try {
        for (std::size_t cell = 0; cell < cells; ++cell) {
          Real pressure_potential = (op.geopotential_from_surface_pressure[level] +
                                     op.reference_specific_volume_m3_kg[level] *
                                         op.pressure_from_surface_pressure[level]) *
                                    perturbation.surface_pressure_pa[cell];
          for (std::size_t source = 0; source < levels; ++source) {
            pressure_potential +=
                op.geopotential_from_potential_temperature_mass[level * levels +
                                                                source] *
                perturbation.potential_temperature_mass_k_kg_m2[dry_hydrostatic_offset(
                    cell, source, levels)];
          }
          worker.geopotential_perturbation[cell] = pressure_potential;
        }
        least_squares_gradient(grid, worker.geopotential_perturbation,
                               worker.geopotential_gradient);
        for (std::size_t cell = 0; cell < cells; ++cell) {
          const auto n = dry_hydrostatic_offset(cell, level, levels);
          result.tendency.momentum[n] =
              -op.reference_air_mass_kg_m2[level] *
              project_tangent(worker.geopotential_gradient[cell],
                              grid.cells()[cell].center);
        }
      } catch (...) {
        workspace.failures[level] = std::current_exception();
      }
    }
    for (std::size_t level = 0; level < levels; ++level)
      if (workspace.failures[level] != nullptr)
        std::rethrow_exception(workspace.failures[level]);
  }
  if (compute_scalar_from_momentum && compute_momentum_from_scalar)
    result.tendency.tracer_mass.assign(volume, 0.0);
}

DryHydrostaticFastTendency subtract_dry_hydrostatic_fast_tendency(
    const DryHydrostaticFastTendency& left, const DryHydrostaticFastTendency& right) {
  return combine_fast_tendencies(left, right, -1.0);
}

DryHydrostaticFastTendency add_dry_hydrostatic_fast_tendency(
    const DryHydrostaticFastTendency& left, const DryHydrostaticFastTendency& right) {
  return combine_fast_tendencies(left, right, 1.0);
}

}  // namespace mps
