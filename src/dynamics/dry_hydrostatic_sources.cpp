#include "myplanetsim/dynamics/dry_hydrostatic_sources.hpp"

#include <stdexcept>

#if defined(MPS_ENABLE_OPENMP)
#include <omp.h>
#endif

#include "myplanetsim/numerics/spherical_operators.hpp"
namespace mps {
namespace {

[[nodiscard]] std::size_t source_worker_count(const std::size_t levels) {
#if defined(MPS_ENABLE_OPENMP)
  return std::min(levels, static_cast<std::size_t>(std::max(1, omp_get_max_threads())));
#else
  static_cast<void>(levels);
  return 1;
#endif
}

[[nodiscard]] std::size_t source_worker_index() {
#if defined(MPS_ENABLE_OPENMP)
  return static_cast<std::size_t>(omp_get_thread_num());
#else
  return 0;
#endif
}

void validate_reference(const CubedSphereGrid& grid, const DryHydrostaticDerived& d,
                        const DryHydrostaticPressureReference& reference) {
  const auto volume = d.cells * d.levels;
  if (reference.cells != grid.cell_count() || reference.cells != d.cells ||
      reference.levels != d.levels || reference.pressure_pa.size() != volume ||
      reference.geopotential_m2_s2.size() != volume ||
      reference.specific_volume_m3_kg.size() != volume ||
      reference.pressure_gradient_pa_m.size() != volume)
    throw std::invalid_argument("dry pressure reference shape mismatch");
}

void calculate_sources(const CubedSphereGrid& grid, const DryHydrostaticDerived& d,
                       const PlanetParameters& p,
                       const DryHydrostaticPressureReference* reference,
                       DryHydrostaticSources& out,
                       DryHydrostaticSourcesWorkspace& workspace) {
  if (d.cells != grid.cell_count() || d.levels == 0)
    throw std::invalid_argument("dry source shape mismatch");
  if (reference != nullptr) validate_reference(grid, d, *reference);
  out.geopotential_gradient_kg_m_s2.resize(d.cells * d.levels);
  out.pressure_correction_kg_m_s2.resize(d.cells * d.levels);
  out.pressure_gradient_kg_m_s2.resize(d.cells * d.levels);
  out.coriolis_kg_m_s2.resize(d.cells * d.levels);
  const Vec3 omega{0, 0, p.rotation_rate_rad_s};
  const auto workers = source_worker_count(d.levels);
  workspace.workers.resize(workers);
  for (auto& worker : workspace.workers) {
    worker.pressure.resize(d.cells);
    worker.geopotential.resize(d.cells);
    worker.pressure_gradient.resize(d.cells);
    worker.geopotential_gradient.resize(d.cells);
  }
  workspace.failures.resize(d.levels);
#if defined(MPS_ENABLE_OPENMP)
#pragma omp parallel for schedule(static) num_threads(workers)
#endif
  for (std::ptrdiff_t level_index = 0;
       level_index < static_cast<std::ptrdiff_t>(d.levels); ++level_index) {
    const auto k = static_cast<std::size_t>(level_index);
    auto& worker = workspace.workers[source_worker_index()];
    workspace.failures[k] = nullptr;
    try {
      for (std::size_t c = 0; c < d.cells; ++c) {
        const auto n = dry_hydrostatic_offset(c, k, d.levels);
        worker.pressure[c] =
            d.pressure_pa[n] - (reference == nullptr ? 0.0 : reference->pressure_pa[n]);
        worker.geopotential[c] =
            d.geopotential_m2_s2[n] -
            (reference == nullptr ? 0.0 : reference->geopotential_m2_s2[n]);
      }
      least_squares_gradient(grid, worker.pressure, worker.pressure_gradient);
      least_squares_gradient(grid, worker.geopotential, worker.geopotential_gradient);
      for (std::size_t c = 0; c < d.cells; ++c) {
        const auto n = dry_hydrostatic_offset(c, k, d.levels);
        const auto centre = grid.cells()[c].center;
        const auto alpha =
            p.gas_constant_j_kg_k * d.temperature_k[n] / d.pressure_pa[n];
        const auto reference_correction =
            reference == nullptr ? Vec3{}
                                 : (alpha - reference->specific_volume_m3_kg[n]) *
                                       reference->pressure_gradient_pa_m[n];
        out.geopotential_gradient_kg_m_s2[n] =
            -d.air_mass_kg_m2[n] *
            project_tangent(worker.geopotential_gradient[c], centre);
        out.pressure_correction_kg_m_s2[n] =
            -d.air_mass_kg_m2[n] *
            project_tangent(alpha * worker.pressure_gradient[c] + reference_correction,
                            centre);
        out.pressure_gradient_kg_m_s2[n] =
            out.geopotential_gradient_kg_m_s2[n] + out.pressure_correction_kg_m_s2[n];
        out.coriolis_kg_m_s2[n] =
            -2 * d.air_mass_kg_m2[n] *
            project_tangent(cross(omega, d.velocity_m_s[n]), centre);
      }
    } catch (...) {
      workspace.failures[k] = std::current_exception();
    }
  }
  for (std::size_t level = 0; level < d.levels; ++level)
    if (workspace.failures[level] != nullptr)
      std::rethrow_exception(workspace.failures[level]);
}

}  // namespace

DryHydrostaticSources dry_hydrostatic_sources(const CubedSphereGrid& grid,
                                              const DryHydrostaticDerived& d,
                                              const PlanetParameters& p) {
  DryHydrostaticSources result;
  DryHydrostaticSourcesWorkspace workspace;
  dry_hydrostatic_sources(grid, d, p, result, workspace);
  return result;
}

DryHydrostaticPressureReference make_dry_hydrostatic_pressure_reference(
    const CubedSphereGrid& grid, const DryHydrostaticDerived& d,
    const PlanetParameters& p) {
  if (d.cells != grid.cell_count() || d.levels == 0 ||
      d.pressure_pa.size() != d.cells * d.levels ||
      d.geopotential_m2_s2.size() != d.cells * d.levels ||
      d.temperature_k.size() != d.cells * d.levels)
    throw std::invalid_argument("dry pressure reference shape mismatch");
  DryHydrostaticPressureReference reference{
      .cells = d.cells,
      .levels = d.levels,
      .pressure_pa = d.pressure_pa,
      .geopotential_m2_s2 = d.geopotential_m2_s2,
      .specific_volume_m3_kg = std::vector<Real>(d.cells * d.levels),
      .pressure_gradient_pa_m = std::vector<Vec3>(d.cells * d.levels),
  };
  std::vector<Real> pressure(d.cells);
  std::vector<Vec3> gradient(d.cells);
  for (std::size_t k = 0; k < d.levels; ++k) {
    for (std::size_t c = 0; c < d.cells; ++c)
      pressure[c] = d.pressure_pa[dry_hydrostatic_offset(c, k, d.levels)];
    least_squares_gradient(grid, pressure, gradient);
    for (std::size_t c = 0; c < d.cells; ++c) {
      const auto n = dry_hydrostatic_offset(c, k, d.levels);
      reference.specific_volume_m3_kg[n] =
          p.gas_constant_j_kg_k * d.temperature_k[n] / d.pressure_pa[n];
      reference.pressure_gradient_pa_m[n] = gradient[c];
    }
  }
  return reference;
}

void dry_hydrostatic_sources(const CubedSphereGrid& grid,
                             const DryHydrostaticDerived& d, const PlanetParameters& p,
                             DryHydrostaticSources& out,
                             DryHydrostaticSourcesWorkspace& workspace) {
  calculate_sources(grid, d, p, nullptr, out, workspace);
}

void dry_hydrostatic_sources(const CubedSphereGrid& grid,
                             const DryHydrostaticDerived& d, const PlanetParameters& p,
                             const DryHydrostaticPressureReference& reference,
                             DryHydrostaticSources& out,
                             DryHydrostaticSourcesWorkspace& workspace) {
  calculate_sources(grid, d, p, &reference, out, workspace);
}
}  // namespace mps
