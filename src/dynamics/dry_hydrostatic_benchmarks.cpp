#include "myplanetsim/dynamics/dry_hydrostatic_benchmarks.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>

namespace mps {
namespace {
constexpr Real kJw06U0 = 35.0;
constexpr Real kJw06Eta0 = 0.252;

Real jw06_vertical_factor(const Real eta) {
  const Real eta_v = (eta - kJw06Eta0) * 0.5 * std::numbers::pi_v<Real>;
  return std::pow(std::max(0.0, std::cos(eta_v)), 1.5);
}

Real jw06_temperature(const Real eta, const Real latitude,
                      const PlanetParameters& planet) {
  constexpr Real reference_temperature = 288.0;
  constexpr Real lapse_rate = 0.005;
  constexpr Real tropopause_eta = 0.2;
  constexpr Real delta_temperature = 4.8e5;
  const Real mean =
      reference_temperature *
          std::pow(eta, planet.gas_constant_j_kg_k * lapse_rate / planet.gravity_m_s2) +
      (eta < tropopause_eta ? delta_temperature * std::pow(tropopause_eta - eta, 5)
                            : 0.0);
  const Real eta_v = (eta - kJw06Eta0) * 0.5 * std::numbers::pi_v<Real>;
  const Real cosine_eta = std::max(0.0, std::cos(eta_v));
  const Real sine = std::sin(latitude);
  const Real cosine = std::cos(latitude);
  const Real wind_shape =
      -2.0 * std::pow(sine, 6) * (cosine * cosine + 1.0 / 3.0) + 10.0 / 63.0;
  const Real rotation_shape =
      8.0 / 5.0 * std::pow(cosine, 3) * (sine * sine + 2.0 / 3.0) -
      std::numbers::pi_v<Real> / 4.0;
  const Real bracket = 2.0 * kJw06U0 * std::pow(cosine_eta, 1.5) * wind_shape +
                       planet.radius_m * planet.rotation_rate_rad_s * rotation_shape;
  return mean + 0.75 * eta * std::numbers::pi_v<Real> * kJw06U0 /
                    planet.gas_constant_j_kg_k * std::sin(eta_v) *
                    std::sqrt(cosine_eta) * bracket;
}
}  // namespace

DryHydrostaticState initialize_dry_hydrostatic_benchmark(
    const ExperimentConfig& config, const CubedSphereGrid& grid,
    const AtmosphericHybridCoordinate& coordinate, const SurfaceOrography& orography) {
  DryHydrostaticState state{.time_s = config.run.start_time_s};
  const auto cells = grid.cell_count();
  const auto levels = coordinate.levels();
  state.surface_pressure_pa.assign(cells, config.vertical.surface_pressure_pa);
  if (orography.surface_geopotential_m2_s2().size() != cells)
    throw std::invalid_argument("benchmark orography shape does not match grid");
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
    if (config.dry_hydrostatic.test_case == DryHydrostaticTestCase::kDcmip200Rest) {
      constexpr Real reference_temperature_k = 300.0;
      constexpr Real lapse_rate_k_m = 0.0065;
      const Real height_m =
          orography.surface_geopotential_m2_s2()[cell] / config.planet.gravity_m_s2;
      state.surface_pressure_pa[cell] =
          config.planet.reference_pressure_pa *
          std::pow(1.0 - lapse_rate_k_m * height_m / reference_temperature_k,
                   config.planet.gravity_m_s2 /
                       (config.planet.gas_constant_j_kg_k * lapse_rate_k_m));
    }
    if (config.dry_hydrostatic.test_case ==
        DryHydrostaticTestCase::kLinearMountainWave) {
      state.surface_pressure_pa[cell] =
          config.planet.reference_pressure_pa *
          std::exp(-orography.surface_geopotential_m2_s2()[cell] /
                   (config.planet.gas_constant_j_kg_k *
                    config.vertical.initial_temperature_k));
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
        case DryHydrostaticTestCase::kJw06Steady:
        case DryHydrostaticTestCase::kJw06Baroclinic: {
          const Real eta =
              geometry.pressure_full_pa[level] / config.planet.reference_pressure_pa;
          Real zonal = kJw06U0 * jw06_vertical_factor(eta) *
                       std::pow(std::sin(2.0 * latitude), 2);
          if (config.dry_hydrostatic.test_case ==
              DryHydrostaticTestCase::kJw06Baroclinic) {
            constexpr Real centre_longitude = std::numbers::pi_v<Real> / 9.0;
            constexpr Real centre_latitude = 2.0 * std::numbers::pi_v<Real> / 9.0;
            const Vec3 centre{std::cos(centre_latitude) * std::cos(centre_longitude),
                              std::cos(centre_latitude) * std::sin(centre_longitude),
                              std::sin(centre_latitude)};
            const Real angular_distance = safe_angle(position, centre);
            zonal += std::exp(-100.0 * angular_distance * angular_distance);
          }
          velocity = zonal * normalize(cross(Vec3{0, 0, 1}, position));
          break;
        }
        case DryHydrostaticTestCase::kLinearMountainWave:
          velocity = 10.0 * cross(Vec3{0, 0, 1}, position);
          break;
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
      Real temperature = config.vertical.initial_temperature_k;
      if (config.dry_hydrostatic.test_case == DryHydrostaticTestCase::kJw06Steady ||
          config.dry_hydrostatic.test_case == DryHydrostaticTestCase::kJw06Baroclinic) {
        const Real eta =
            geometry.pressure_full_pa[level] / config.planet.reference_pressure_pa;
        temperature = jw06_temperature(eta, latitude, config.planet);
      }
      if (config.dry_hydrostatic.test_case == DryHydrostaticTestCase::kDcmip200Rest) {
        constexpr Real lapse_rate_k_m = 0.0065;
        temperature = 300.0 * std::pow(geometry.pressure_full_pa[level] /
                                           config.planet.reference_pressure_pa,
                                       config.planet.gas_constant_j_kg_k *
                                           lapse_rate_k_m / config.planet.gravity_m_s2);
      }
      const Real theta = temperature / geometry.exner_full[level];
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
