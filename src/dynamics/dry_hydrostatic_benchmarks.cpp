#include "myplanetsim/dynamics/dry_hydrostatic_benchmarks.hpp"

#include <cmath>

namespace mps {

DryHydrostaticState initialize_dry_hydrostatic_benchmark(
    const ExperimentConfig& config, const CubedSphereGrid& grid,
    const AtmosphericHybridCoordinate& coordinate) {
  DryHydrostaticState state{.time_s = config.run.start_time_s};
  const auto cells = grid.cell_count();
  const auto levels = coordinate.levels();
  state.surface_pressure_pa.assign(cells, config.vertical.surface_pressure_pa);
  state.horizontal_momentum_mass_kg_m_s.resize(cells * levels);
  state.potential_temperature_mass_k_kg_m2.resize(cells * levels);
  state.tracer_mass_kg_m2.resize(cells * levels);

  for (std::size_t cell = 0; cell < cells; ++cell) {
    const auto position = grid.cells()[cell].center;
    const Real longitude = std::atan2(position.y, position.x);
    const Real latitude = std::asin(position.z);
    if (config.dry_hydrostatic.test_case == DryHydrostaticTestCase::kLinearWave) {
      state.surface_pressure_pa[cell] *=
          1.0 + 1.0e-5 * std::cos(longitude) * std::cos(latitude);
    }
    const auto geometry = coordinate.geometry(
        state.surface_pressure_pa[cell], config.planet.gravity_m_s2,
        config.planet.gas_constant_j_kg_k, config.planet.heat_capacity_cp_j_kg_k,
        config.planet.reference_pressure_pa);
    for (std::size_t level = 0; level < levels; ++level) {
      const auto offset = dry_hydrostatic_offset(cell, level, levels);
      const Real sigma = (static_cast<Real>(level) + 0.5) / static_cast<Real>(levels);
      Vec3 velocity{};
      Real tracer = 0.0;
      switch (config.dry_hydrostatic.test_case) {
        case DryHydrostaticTestCase::kSolidBodyTransport:
          velocity = 20.0 * cross(Vec3{0, 0, 1}, position);
          tracer = std::exp(-20.0 * ((longitude - 0.5) * (longitude - 0.5) +
                                     latitude * latitude)) *
                   (1.0 - sigma);
          break;
        case DryHydrostaticTestCase::kDcmipDeformational:
          velocity = 10.0 * std::sin(2.0 * longitude) *
                     project_tangent(Vec3{0, 0, 1}, position);
          tracer = 0.5 + 0.25 * std::cos(longitude) * std::cos(latitude) *
                             std::sin(3.14159265358979323846 * sigma);
          break;
        case DryHydrostaticTestCase::kDcmipHadley:
          velocity =
              8.0 * std::sin(2.0 * latitude) * project_tangent(Vec3{0, 0, 1}, position);
          tracer = 0.5 + 0.2 * position.z * (1.0 - sigma);
          break;
        case DryHydrostaticTestCase::kUmjs14Steady:
        case DryHydrostaticTestCase::kUmjs14Baroclinic:
          velocity =
              35.0 * std::cos(latitude) * normalize(cross(Vec3{0, 0, 1}, position));
          break;
        default:
          break;
      }
      const Real theta =
          config.vertical.initial_temperature_k / geometry.exner_full[level];
      state.horizontal_momentum_mass_kg_m_s[offset] =
          geometry.air_mass_kg_m2[level] * velocity;
      state.potential_temperature_mass_k_kg_m2[offset] =
          geometry.air_mass_kg_m2[level] * theta;
      state.tracer_mass_kg_m2[offset] = geometry.air_mass_kg_m2[level] * tracer;
    }
  }
  return state;
}

}  // namespace mps
