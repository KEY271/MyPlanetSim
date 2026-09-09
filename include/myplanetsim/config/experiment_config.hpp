#pragma once

#include <cstdint>
#include <filesystem>
#include <iosfwd>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include "myplanetsim/core/orbit_parameters.hpp"
#include "myplanetsim/core/planet_parameters.hpp"
#include "myplanetsim/core/types.hpp"
#include "myplanetsim/dynamics/tracer_registry.hpp"

namespace mps {

struct RunParameters {
  Real start_time_s = 0.0;
  Real end_time_s = 0.0;
  Real time_step_s = 0.0;
  Seed random_seed = 0;
};

struct OdeParameters {
  Real initial_value = 0.0;
  Real decay_rate_s_1 = 0.0;
};

enum class ExperimentKind { kOde, kSphereTransport, kVerticalColumn, kDryHydrostatic };
enum class TransportScheme { kUpwind, kLinear };
enum class LimiterKind { kNone, kBarthJespersen };
enum class TransportTestCase { kSolidBody, kDeformational, kDivergent };
enum class ReconstructionKind { kPiecewiseConstant, kLinear };
enum class DiffusionKind { kNone, kLaplacian, kBiharmonic };
enum class PhysicsKind {
  kNone,
  kHeldSuarez,
  kPlanetaryNewtonian,
  kSurfaceEnergyBalance,
  kGrayRadiation
};
enum class ConvectionKind { kNone, kDryAdjustment, kSimpleBettsMiller };
enum class BoundaryLayerKind { kNone, kBulkKProfile };
enum class BoundaryLayerIntegrator { kBackwardEuler };
enum class MoistureKind { kNone, kDiluteWater };
enum class CondensationKind { kNone, kSaturationAdjustment };
enum class SurfaceMoistureExchange { kNone, kBulk };
enum class SurfaceHydrologyKind { kNone, kBucket };
enum class ForcingGeometry { kAxisymmetric, kSubstellar };
enum class OrographyKind { kFlat, kDcmip200, kLinearBell, kJw06, kLatLonCsv };
enum class SurfaceGeography { kUniform, kEarth };
enum class InitialConditionKind {
  kConstant,
  kGaussianHill,
  kCosineBell,
  kSlottedCylinder
};
enum class VerticalTestCase {
  kIsothermal,
  kDryAdiabatic,
  kMovingSurfacePressure,
  kManufacturedTransport
};
enum class VerticalTransportScheme { kDonorCell, kLinear };
enum class VerticalLimiterKind { kNone, kMinmod };
enum class DryHydrostaticTimeIntegrator {
  kExplicitSspRk3,
  kSemiImplicit,
  kArk2ImexComparison
};
enum class DryHydrostaticTestCase {
  kIsothermalRest,
  kSolidBodyTransport,
  kDcmipDeformational,
  kDcmipHadley,
  kLinearWave,
  kDcmip200Rest,
  kLinearMountainWave,
  kJw06Steady,
  kJw06Baroclinic,
  kUmjs14Steady,
  kUmjs14Baroclinic,
  kHeldSuarez
};

struct GridParameters {
  Index cells_per_panel = 0;
};

struct TransportParameters {
  TransportTestCase test_case = TransportTestCase::kSolidBody;
  InitialConditionKind initial_condition = InitialConditionKind::kGaussianHill;
  TransportScheme scheme = TransportScheme::kLinear;
  LimiterKind limiter = LimiterKind::kBarthJespersen;
  Real cfl = 0.5;
  Real rotation_axis_x = 0.0;
  Real rotation_axis_y = 0.0;
  Real rotation_axis_z = 1.0;
  Real angular_speed_rad_s = 0.0;
};

struct VerticalColumnParameters {
  VerticalTestCase test_case = VerticalTestCase::kIsothermal;
  Index levels = 0;
  std::vector<Real> a_half_pa{};
  std::vector<Real> b_half{};
  Real surface_pressure_pa = 0.0;
  Real minimum_surface_pressure_pa = 0.0;
  Real maximum_surface_pressure_pa = 0.0;
  Real minimum_pressure_thickness_pa = 0.0;
  Real surface_geopotential_m2_s2 = 0.0;
  Real initial_temperature_k = 0.0;
  Real initial_potential_temperature_k = 0.0;
  Real temperature_floor_k = 0.0;
  VerticalTransportScheme transport_scheme = VerticalTransportScheme::kDonorCell;
  VerticalLimiterKind limiter = VerticalLimiterKind::kNone;
  Real cfl = 0.5;
  Real forcing_amplitude = 0.0;
};

struct DryHydrostaticParameters {
  DryHydrostaticTestCase test_case = DryHydrostaticTestCase::kIsothermalRest;
  ReconstructionKind reconstruction = ReconstructionKind::kLinear;
  LimiterKind limiter = LimiterKind::kBarthJespersen;
  Real cfl = 0.45;
  DiffusionKind diffusion_kind = DiffusionKind::kNone;
  Real diffusion_coefficient = 0.0;
  DryHydrostaticTimeIntegrator time_integrator =
      DryHydrostaticTimeIntegrator::kExplicitSspRk3;
  Real advective_cfl = 0.45;
};

// ADR 0016: the reference-linear operator may be frozen at construction or rebuilt
// from the horizontal mean of the current state at every step. The per-step form
// follows Thuburn et al. (2014), whose reference field is updated from the state at
// step n and whose quasi-Newton residual then falls about one order per iteration.
enum class SemiImplicitReferenceUpdate { kFixed, kPerStep };

struct SemiImplicitParameters {
  Real reference_surface_pressure_pa = 0.0;
  Real reference_temperature_k = 0.0;
  SemiImplicitReferenceUpdate reference_update = SemiImplicitReferenceUpdate::kFixed;
  Real implicit_weight = 0.5;
  Real wave_cfl_threshold = 0.45;
  Index maximum_implicit_modes = 0;
  // ICI schemes are indexed by a fixed iteration count, not by a residual tolerance
  // (Benard 2003). The residual is recorded as a diagnostic only.
  Index nonlinear_iterations = 0;
  Real linear_relative_tolerance = 0.0;
  Real linear_absolute_tolerance = 0.0;
  Index linear_maximum_iterations = 0;
  Index gmres_restart = 0;
  Real minimum_time_step_s = 0.0;
};

struct OrographyParameters {
  OrographyKind kind = OrographyKind::kFlat;
  std::string input_file{};
  std::string input_fingerprint_fnv1a64{};
  Index smoothing_passes = 0;
};

struct DiagnosticsParameters {
  std::uint64_t interval_steps = 1;
};

struct StatisticsParameters {
  bool enabled = false;
  Real start_time_s = 0.0;
  Real period_s = 0.0;
  std::vector<std::string> fields{};
};

struct PhysicsParameters {
  PhysicsKind kind = PhysicsKind::kNone;
  ForcingGeometry geometry = ForcingGeometry::kAxisymmetric;
};

// Legacy keeps the Phase 13 common substep and arithmetic order. Process intervals
// enables the Phase 14 event scheduler; its intervals are maxima and are shortened at
// every accepted dynamics-step boundary.
enum class PhysicsScheduleKind { kLegacy, kProcessIntervals };
enum class ConvectionUpdateMode { kIntermittent, kCachedRelaxation };

struct PhysicsScheduleParameters {
  PhysicsScheduleKind kind = PhysicsScheduleKind::kLegacy;
  Real boundary_layer_maximum_update_interval_s = 300.0;
  Real convection_diagnostic_interval_s = 300.0;
  ConvectionUpdateMode convection_update_mode = ConvectionUpdateMode::kIntermittent;
  Real radiation_diagnostic_interval_s = 300.0;
};

struct ConvectionParameters {
  ConvectionKind kind = ConvectionKind::kNone;
  Real stability_tolerance_k = 1e-10;
  Real relaxation_time_s = 7200.0;
  Real reference_relative_humidity = 0.8;
};

struct MoistureParameters {
  MoistureKind kind = MoistureKind::kNone;
  CondensationKind condensation = CondensationKind::kNone;
  SurfaceMoistureExchange surface_exchange = SurfaceMoistureExchange::kNone;
  Real maximum_physics_substep_s = 300.0;
};

struct BoundaryLayerParameters {
  BoundaryLayerKind kind = BoundaryLayerKind::kNone;
  BoundaryLayerIntegrator integrator = BoundaryLayerIntegrator::kBackwardEuler;
  Real critical_richardson = 1.0;
  Real turbulent_prandtl = 1.0;
  Real gustiness_m_s = 1.0;
};

struct SurfaceParameters {
  SurfaceGeography geography = SurfaceGeography::kUniform;
  Real uniform_land_fraction = 0.0;
  std::string input_file{};
  std::string input_fingerprint_fnv1a64{};
  Index quadrature_order = 1;
  Index smoothing_passes = 0;
  Real land_heat_capacity_j_m2_k = 0.0;
  Real ocean_heat_capacity_j_m2_k = 0.0;
  Real initial_temperature_k = 0.0;
  Real albedo = 0.0;
  Real emissivity = 0.0;
  Real air_exchange_coefficient_w_m2_k = 0.0;
  Real internal_heat_flux_w_m2 = 0.0;
  Real land_roughness_momentum_m = 0.0;
  Real land_roughness_heat_m = 0.0;
  Real ocean_roughness_momentum_m = 0.0;
  Real ocean_roughness_heat_m = 0.0;
  SurfaceHydrologyKind hydrology_kind = SurfaceHydrologyKind::kNone;
  Real hydrology_capacity_kg_m2 = 0.0;
  Real hydrology_initial_fraction = 0.0;
  // Stability fraction for the explicit surface reservoir (ADR 0011):
  // dt <= cfl * min_c C_surface[c] / (4 eps sigma_SB T_s[c]^3).
  Real cfl = 0.5;
};

struct RadiationParameters {
  Real shortwave_absorption_m2_kg = 0.0;
  Real longwave_absorption_ref_m2_kg = 0.0;
  Real reference_pressure_pa = 0.0;
  Real longwave_pressure_exponent = 1.0;
  Real longwave_diffusivity_factor = 1.66;
  Real shortwave_diffuse_factor = 1.66;
  Real cfl = 0.5;
};

struct ExperimentConfig {
  ExperimentKind kind = ExperimentKind::kOde;
  PlanetParameters planet{};
  RunParameters run{};
  OdeParameters ode{};
  GridParameters grid{};
  TransportParameters transport{};
  VerticalColumnParameters vertical{};
  DryHydrostaticParameters dry_hydrostatic{};
  std::vector<TracerDescriptor> tracers{};
  std::optional<SemiImplicitParameters> semi_implicit{};
  OrographyParameters orography{};
  PhysicsParameters physics{};
  PhysicsScheduleParameters physics_schedule{};
  ConvectionParameters convection{};
  MoistureParameters moisture{};
  BoundaryLayerParameters boundary_layer{};
  DiagnosticsParameters diagnostics{};
  StatisticsParameters statistics{};
  std::string output_directory{};
  // Runtime-only origin used to resolve portable config-relative inputs.
  std::filesystem::path source_directory{};
  std::optional<OrbitParameters> orbit{};
  std::optional<SurfaceParameters> surface{};
  std::optional<RadiationParameters> radiation{};

  void validate() const;
};

[[nodiscard]] std::string_view experiment_kind_name(ExperimentKind kind) noexcept;
[[nodiscard]] std::string_view transport_scheme_name(TransportScheme scheme) noexcept;
[[nodiscard]] std::string_view limiter_name(LimiterKind limiter) noexcept;
[[nodiscard]] std::string_view transport_test_case_name(
    TransportTestCase test_case) noexcept;
[[nodiscard]] std::string_view initial_condition_name(
    InitialConditionKind initial_condition) noexcept;
[[nodiscard]] std::string_view reconstruction_name(
    ReconstructionKind reconstruction) noexcept;
[[nodiscard]] std::string_view diffusion_kind_name(DiffusionKind kind) noexcept;
[[nodiscard]] std::string_view vertical_test_case_name(
    VerticalTestCase test_case) noexcept;
[[nodiscard]] std::string_view vertical_transport_scheme_name(
    VerticalTransportScheme scheme) noexcept;
[[nodiscard]] std::string_view vertical_limiter_name(
    VerticalLimiterKind limiter) noexcept;
[[nodiscard]] std::string_view dry_hydrostatic_test_case_name(
    DryHydrostaticTestCase test_case) noexcept;
[[nodiscard]] std::string_view dry_hydrostatic_time_integrator_name(
    DryHydrostaticTimeIntegrator integrator) noexcept;
[[nodiscard]] std::string_view orography_kind_name(OrographyKind kind) noexcept;
[[nodiscard]] std::string_view physics_kind_name(PhysicsKind kind) noexcept;
[[nodiscard]] std::string_view physics_schedule_kind_name(
    PhysicsScheduleKind kind) noexcept;
[[nodiscard]] std::string_view convection_update_mode_name(
    ConvectionUpdateMode mode) noexcept;
[[nodiscard]] std::string_view convection_kind_name(ConvectionKind kind) noexcept;
[[nodiscard]] std::string_view moisture_kind_name(MoistureKind kind) noexcept;
[[nodiscard]] std::string_view boundary_layer_kind_name(
    BoundaryLayerKind kind) noexcept;
[[nodiscard]] std::string_view boundary_layer_integrator_name(
    BoundaryLayerIntegrator integrator) noexcept;
[[nodiscard]] std::string_view forcing_geometry_name(ForcingGeometry geometry) noexcept;
[[nodiscard]] std::string_view surface_geography_name(
    SurfaceGeography geography) noexcept;

[[nodiscard]] ExperimentConfig parse_experiment_config(std::istream& input);
[[nodiscard]] ExperimentConfig load_experiment_config(
    const std::filesystem::path& path);
void write_experiment_config(std::ostream& output, const ExperimentConfig& config);

}  // namespace mps
