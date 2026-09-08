#include "myplanetsim/diagnostics/moist_diagnostics.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace mps {
namespace {

void require_shapes(const CubedSphereGrid& grid, const DryHydrostaticState& state,
                    const DryHydrostaticDerived& derived,
                    const TracerRegistry& registry) {
  const std::size_t volume = derived.cells * derived.levels;
  if (derived.cells != grid.cell_count() || derived.levels == 0 ||
      state.tracer_count != registry.size() ||
      derived.tracer_count != registry.size() ||
      state.tracer_mass_kg_m2.size() != registry.size() * volume ||
      derived.tracer_mixing_ratio.size() != registry.size() * volume ||
      derived.air_mass_kg_m2.size() != volume || derived.pressure_pa.size() != volume ||
      derived.temperature_k.size() != volume)
    throw std::invalid_argument("moist diagnostics shape mismatch");
}

}  // namespace

std::vector<TracerGlobalDiagnostics> diagnose_tracers(
    const CubedSphereGrid& grid, const DryHydrostaticState& state,
    const DryHydrostaticDerived& derived, const TracerRegistry& registry) {
  require_shapes(grid, state, derived, registry);
  std::vector<TracerGlobalDiagnostics> result(registry.size());
  for (std::size_t tracer = 0; tracer < registry.size(); ++tracer) {
    auto& diagnostic = result[tracer];
    diagnostic.minimum_mixing_ratio = std::numeric_limits<Real>::infinity();
    diagnostic.maximum_mixing_ratio = -std::numeric_limits<Real>::infinity();
    for (std::size_t cell = 0; cell < derived.cells; ++cell) {
      const Real area = grid.cells()[cell].area_m2;
      for (std::size_t level = 0; level < derived.levels; ++level) {
        const std::size_t offset = dry_hydrostatic_tracer_offset(
            tracer, cell, level, derived.cells, derived.levels);
        diagnostic.mass_kg += area * state.tracer_mass_kg_m2[offset];
        diagnostic.minimum_mixing_ratio = std::min(diagnostic.minimum_mixing_ratio,
                                                   derived.tracer_mixing_ratio[offset]);
        diagnostic.maximum_mixing_ratio = std::max(diagnostic.maximum_mixing_ratio,
                                                   derived.tracer_mixing_ratio[offset]);
      }
    }
  }
  return result;
}

MoistureGlobalDiagnostics diagnose_moisture(
    const CubedSphereGrid& grid, const DryHydrostaticState& state,
    const DryHydrostaticDerived& derived, const std::span<const Real> land_fraction,
    const std::size_t water_vapor_tracer,
    const DiluteMoistThermodynamics& thermodynamics,
    const Real dilute_limit_mixing_ratio) {
  const std::size_t volume = derived.cells * derived.levels;
  if (derived.cells != grid.cell_count() || derived.levels == 0 ||
      water_vapor_tracer >= state.tracer_count ||
      state.tracer_count != derived.tracer_count ||
      state.tracer_mass_kg_m2.size() != state.tracer_count * volume ||
      derived.tracer_mixing_ratio.size() != derived.tracer_count * volume ||
      derived.air_mass_kg_m2.size() != volume || derived.pressure_pa.size() != volume ||
      derived.temperature_k.size() != volume ||
      state.land_water_kg_m2.size() != derived.cells ||
      land_fraction.size() != derived.cells || !(dilute_limit_mixing_ratio > 0.0))
    throw std::invalid_argument("moist diagnostics shape or limit is invalid");
  MoistureGlobalDiagnostics result;
  result.minimum_relative_humidity = std::numeric_limits<Real>::infinity();
  Real relative_humidity_area_mass_sum = 0.0;
  Real air_area_mass_sum = 0.0;
  Real dilute_exceedance_area = 0.0;
  for (std::size_t cell = 0; cell < derived.cells; ++cell) {
    const Real area = grid.cells()[cell].area_m2;
    if (!std::isfinite(land_fraction[cell]) || land_fraction[cell] < 0.0 ||
        land_fraction[cell] > 1.0)
      throw std::invalid_argument("moist diagnostics land fraction is invalid");
    if (!std::isfinite(state.land_water_kg_m2[cell]) ||
        state.land_water_kg_m2[cell] < 0.0)
      throw std::invalid_argument("moist diagnostics land water is invalid");
    bool cell_exceeds_dilute_limit = false;
    result.land_water_kg += area * land_fraction[cell] * state.land_water_kg_m2[cell];
    for (std::size_t level = 0; level < derived.levels; ++level) {
      const std::size_t scalar = dry_hydrostatic_offset(cell, level, derived.levels);
      const std::size_t vapor = dry_hydrostatic_tracer_offset(
          water_vapor_tracer, cell, level, derived.cells, derived.levels);
      const Real mixing_ratio = derived.tracer_mixing_ratio[vapor];
      const Real relative =
          relative_humidity(mixing_ratio, derived.temperature_k[scalar],
                            derived.pressure_pa[scalar], thermodynamics);
      const Real air_mass = area * derived.air_mass_kg_m2[scalar];
      result.atmospheric_water_kg += area * state.tracer_mass_kg_m2[vapor];
      relative_humidity_area_mass_sum += air_mass * relative;
      air_area_mass_sum += air_mass;
      result.minimum_relative_humidity =
          std::min(result.minimum_relative_humidity, relative);
      result.maximum_relative_humidity =
          std::max(result.maximum_relative_humidity, relative);
      result.maximum_vapor_mixing_ratio =
          std::max(result.maximum_vapor_mixing_ratio, mixing_ratio);
      cell_exceeds_dilute_limit =
          cell_exceeds_dilute_limit || mixing_ratio > dilute_limit_mixing_ratio;
    }
    if (cell_exceeds_dilute_limit) dilute_exceedance_area += area;
  }
  const Real total_area = grid.total_area_m2();
  if (!(total_area > 0.0) || !(air_area_mass_sum > 0.0))
    throw std::runtime_error("moist diagnostics normalization is invalid");
  result.precipitable_water_kg_m2 = result.atmospheric_water_kg / total_area;
  result.area_mean_relative_humidity =
      relative_humidity_area_mass_sum / air_area_mass_sum;
  result.area_mean_land_water_kg_m2 = result.land_water_kg / total_area;
  result.dilute_limit_exceedance_area_fraction = dilute_exceedance_area / total_area;
  return result;
}

}  // namespace mps
