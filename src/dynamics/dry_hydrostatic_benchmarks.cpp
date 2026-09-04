#include "myplanetsim/dynamics/dry_hydrostatic_benchmarks.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <numbers>
#include <stdexcept>
#include <vector>

#include "myplanetsim/dynamics/jw06_parameters.hpp"

namespace mps {
namespace {
[[nodiscard]] std::uint64_t splitmix64(std::uint64_t value) noexcept {
  value += 0x9e3779b97f4a7c15ULL;
  value = (value ^ (value >> 30U)) * 0xbf58476d1ce4e5b9ULL;
  value = (value ^ (value >> 27U)) * 0x94d049bb133111ebULL;
  return value ^ (value >> 31U);
}

[[nodiscard]] Real held_suarez_perturbation(const Seed seed, const std::size_t cell,
                                            const std::size_t level) noexcept {
  const auto key = static_cast<std::uint64_t>(seed) ^
                   (static_cast<std::uint64_t>(cell) * 0xd2b74407b1ce6e93ULL) ^
                   (static_cast<std::uint64_t>(level) * 0xca5a826395121157ULL);
  const Real uniform = static_cast<Real>(splitmix64(key) >> 11U) * 0x1.0p-53;
  return 0.1 * (uniform - 0.5);
}

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
  const Real bracket = 2.0 * kJw06U0Mps * std::pow(cosine_eta, 1.5) * wind_shape +
                       planet.radius_m * planet.rotation_rate_rad_s * rotation_shape;
  return mean + 0.75 * eta * std::numbers::pi_v<Real> * kJw06U0Mps /
                    planet.gas_constant_j_kg_k * std::sin(eta_v) *
                    std::sqrt(cosine_eta) * bracket;
}

struct Umjs14BaseState {
  Real pressure_pa;
  Real temperature_k;
};

[[nodiscard]] Umjs14BaseState umjs14_pressure_temperature(
    const Real height_m, const Real latitude, const PlanetParameters& planet) {
  constexpr Real equatorial_temperature_k = 310.0;
  constexpr Real polar_temperature_k = 240.0;
  constexpr Real jet_half_width = 2.0;
  constexpr Real jet_width = 3.0;
  constexpr Real lapse_rate_k_m = 0.005;
  constexpr Real mean_temperature_k =
      0.5 * (equatorial_temperature_k + polar_temperature_k);
  const Real scaled_height =
      height_m * planet.gravity_m_s2 /
      (jet_half_width * planet.gas_constant_j_kg_k * mean_temperature_k);
  const Real gaussian = std::exp(-scaled_height * scaled_height);
  const Real tau1 =
      std::exp(lapse_rate_k_m * height_m / mean_temperature_k) / mean_temperature_k +
      (mean_temperature_k - polar_temperature_k) /
          (mean_temperature_k * polar_temperature_k) *
          (1.0 - 2.0 * scaled_height * scaled_height) * gaussian;
  const Real tau2 = 0.5 * (jet_width + 2.0) *
                    (equatorial_temperature_k - polar_temperature_k) /
                    (equatorial_temperature_k * polar_temperature_k) *
                    (1.0 - 2.0 * scaled_height * scaled_height) * gaussian;
  const Real cosine = std::cos(latitude);
  const Real temperature_shape =
      std::pow(cosine, jet_width) -
      jet_width / (jet_width + 2.0) * std::pow(cosine, jet_width + 2.0);
  const Real temperature_k = 1.0 / (tau1 - tau2 * temperature_shape);
  const Real integrated_tau1 =
      (std::exp(lapse_rate_k_m * height_m / mean_temperature_k) - 1.0) /
          lapse_rate_k_m +
      height_m * (mean_temperature_k - polar_temperature_k) /
          (mean_temperature_k * polar_temperature_k) * gaussian;
  const Real integrated_tau2 =
      0.5 * (jet_width + 2.0) * (equatorial_temperature_k - polar_temperature_k) /
      (equatorial_temperature_k * polar_temperature_k) * height_m * gaussian;
  const Real pressure_pa =
      planet.reference_pressure_pa *
      std::exp(-planet.gravity_m_s2 / planet.gas_constant_j_kg_k *
               (integrated_tau1 - integrated_tau2 * temperature_shape));
  return {pressure_pa, temperature_k};
}

[[nodiscard]] Real umjs14_height_for_pressure(const Real pressure_pa,
                                              const Real latitude,
                                              const PlanetParameters& planet) {
  Real lower_height = 0.0;
  Real upper_height = 10000.0;
  Real lower_pressure =
      umjs14_pressure_temperature(lower_height, latitude, planet).pressure_pa;
  Real upper_pressure =
      umjs14_pressure_temperature(upper_height, latitude, planet).pressure_pa;
  Real height = upper_height;
  for (int iteration = 0; iteration < 100; ++iteration) {
    const Real denominator = upper_pressure - lower_pressure;
    if (denominator == 0.0)
      throw std::runtime_error("UMJS14 pressure inversion did not converge");
    height = upper_height - (upper_pressure - pressure_pa) *
                                (upper_height - lower_height) / denominator;
    const Real candidate_pressure =
        umjs14_pressure_temperature(height, latitude, planet).pressure_pa;
    if (std::abs((candidate_pressure - pressure_pa) / pressure_pa) < 1.0e-13)
      return height;
    lower_height = upper_height;
    lower_pressure = upper_pressure;
    upper_height = height;
    upper_pressure = candidate_pressure;
  }
  throw std::runtime_error("UMJS14 pressure inversion did not converge");
}

[[nodiscard]] Real umjs14_zonal_wind(const Real height_m, const Real latitude,
                                     const Real temperature_k,
                                     const PlanetParameters& planet) {
  constexpr Real equatorial_temperature_k = 310.0;
  constexpr Real polar_temperature_k = 240.0;
  constexpr Real jet_half_width = 2.0;
  constexpr Real jet_width = 3.0;
  constexpr Real mean_temperature_k =
      0.5 * (equatorial_temperature_k + polar_temperature_k);
  const Real scaled_height =
      height_m * planet.gravity_m_s2 /
      (jet_half_width * planet.gas_constant_j_kg_k * mean_temperature_k);
  const Real integrated_tau2 = 0.5 * (jet_width + 2.0) *
                               (equatorial_temperature_k - polar_temperature_k) /
                               (equatorial_temperature_k * polar_temperature_k) *
                               height_m * std::exp(-scaled_height * scaled_height);
  const Real cosine = std::cos(latitude);
  const Real wind_shape =
      std::pow(cosine, jet_width - 1.0) - std::pow(cosine, jet_width + 1.0);
  const Real forcing = planet.gravity_m_s2 / planet.radius_m * jet_width *
                       integrated_tau2 * wind_shape * temperature_k;
  const Real rotation_speed = planet.rotation_rate_rad_s * planet.radius_m * cosine;
  return -rotation_speed + std::sqrt(rotation_speed * rotation_speed +
                                     planet.radius_m * cosine * forcing);
}

[[nodiscard]] Real umjs14_streamfunction(const Real longitude, const Real latitude,
                                         const Real height_m) {
  constexpr Real perturbation_speed_m_s = 0.5;
  constexpr Real perturbation_radius = 1.0 / 6.0;
  constexpr Real centre_longitude = std::numbers::pi_v<Real> / 9.0;
  constexpr Real centre_latitude = 2.0 * std::numbers::pi_v<Real> / 9.0;
  constexpr Real height_cap_m = 15000.0;
  const Vec3 position{std::cos(latitude) * std::cos(longitude),
                      std::cos(latitude) * std::sin(longitude), std::sin(latitude)};
  const Vec3 centre{std::cos(centre_latitude) * std::cos(centre_longitude),
                    std::cos(centre_latitude) * std::sin(centre_longitude),
                    std::sin(centre_latitude)};
  const Real radius = safe_angle(position, centre) / perturbation_radius;
  if (!(radius < 1.0) || !(height_m < height_cap_m)) return 0.0;
  const Real height_ratio = height_m / height_cap_m;
  const Real vertical_taper = 1.0 - 3.0 * height_ratio * height_ratio +
                              2.0 * height_ratio * height_ratio * height_ratio;
  const Real horizontal_taper = std::cos(0.5 * std::numbers::pi_v<Real> * radius);
  const Real taper_squared = horizontal_taper * horizontal_taper;
  return -perturbation_speed_m_s * perturbation_radius * vertical_taper *
         taper_squared * taper_squared;
}

[[nodiscard]] Vec3 umjs14_velocity(const Real longitude, const Real latitude,
                                   const Real height_m, const Real temperature_k,
                                   const PlanetParameters& planet,
                                   const bool perturbed) {
  Real zonal = umjs14_zonal_wind(height_m, latitude, temperature_k, planet);
  Real meridional = 0.0;
  if (perturbed) {
    constexpr Real epsilon = 1.0e-5;
    zonal -= (umjs14_streamfunction(longitude, latitude + epsilon, height_m) -
              umjs14_streamfunction(longitude, latitude - epsilon, height_m)) /
             (2.0 * epsilon);
    meridional += (umjs14_streamfunction(longitude + epsilon, latitude, height_m) -
                   umjs14_streamfunction(longitude - epsilon, latitude, height_m)) /
                  (2.0 * epsilon * std::cos(latitude));
  }
  const Vec3 position{std::cos(latitude) * std::cos(longitude),
                      std::cos(latitude) * std::sin(longitude), std::sin(latitude)};
  const Vec3 zonal_direction = normalize(cross(Vec3{0, 0, 1}, position));
  const Vec3 north_direction = normalize(project_tangent(Vec3{0, 0, 1}, position));
  return zonal * zonal_direction + meridional * north_direction;
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

  std::vector<Real> held_suarez_temperature_perturbation(cells * levels, 0.0);
  if (config.dry_hydrostatic.test_case == DryHydrostaticTestCase::kHeldSuarez) {
    for (std::size_t level = 0; level < levels; ++level) {
      Real weighted_sum = 0.0;
      Real total_area = 0.0;
      for (std::size_t cell = 0; cell < cells; ++cell) {
        const auto offset = dry_hydrostatic_offset(cell, level, levels);
        held_suarez_temperature_perturbation[offset] =
            held_suarez_perturbation(config.run.random_seed, cell, level);
        weighted_sum +=
            grid.cells()[cell].area_m2 * held_suarez_temperature_perturbation[offset];
        total_area += grid.cells()[cell].area_m2;
      }
      const Real mean = weighted_sum / total_area;
      for (std::size_t cell = 0; cell < cells; ++cell) {
        held_suarez_temperature_perturbation[dry_hydrostatic_offset(cell, level,
                                                                    levels)] -= mean;
      }
    }
  }

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
    if (config.dry_hydrostatic.test_case == DryHydrostaticTestCase::kUmjs14Steady ||
        config.dry_hydrostatic.test_case == DryHydrostaticTestCase::kUmjs14Baroclinic) {
      state.surface_pressure_pa[cell] = config.planet.reference_pressure_pa;
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
          Real zonal = kJw06U0Mps * jw06_vertical_factor(eta) *
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
        case DryHydrostaticTestCase::kUmjs14Baroclinic: {
          const Real height_m = umjs14_height_for_pressure(
              geometry.pressure_full_pa[level], latitude, config.planet);
          const Real balanced_temperature =
              umjs14_pressure_temperature(height_m, latitude, config.planet)
                  .temperature_k;
          velocity = umjs14_velocity(longitude, latitude, height_m,
                                     balanced_temperature, config.planet,
                                     config.dry_hydrostatic.test_case ==
                                         DryHydrostaticTestCase::kUmjs14Baroclinic);
          break;
        }
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
      if (config.dry_hydrostatic.test_case == DryHydrostaticTestCase::kUmjs14Steady ||
          config.dry_hydrostatic.test_case ==
              DryHydrostaticTestCase::kUmjs14Baroclinic) {
        const Real height_m = umjs14_height_for_pressure(
            geometry.pressure_full_pa[level], latitude, config.planet);
        temperature = umjs14_pressure_temperature(height_m, latitude, config.planet)
                          .temperature_k;
      }
      if (config.dry_hydrostatic.test_case == DryHydrostaticTestCase::kHeldSuarez) {
        temperature = 264.0 + held_suarez_temperature_perturbation[offset];
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
