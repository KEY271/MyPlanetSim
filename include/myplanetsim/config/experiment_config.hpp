#pragma once

#include <cstdint>
#include <filesystem>
#include <iosfwd>
#include <string>
#include <string_view>
#include <vector>

#include "myplanetsim/core/planet_parameters.hpp"
#include "myplanetsim/core/types.hpp"

namespace mps {

struct RunParameters {
  Real start_time_s;
  Real end_time_s;
  Real time_step_s;
  Seed random_seed;
};

struct OdeParameters {
  Real initial_value;
  Real decay_rate_s_1;
};

enum class ExperimentKind { kOde, kSphereTransport, kShallowWater, kVerticalColumn };
enum class TransportScheme { kUpwind, kLinear };
enum class LimiterKind { kNone, kBarthJespersen };
enum class TransportTestCase { kSolidBody, kDeformational, kDivergent };
enum class ShallowWaterScheme { kRusanov, kCompatible };
enum class ReconstructionKind { kPiecewiseConstant, kLinear };
enum class ShallowWaterTestCase {
  kRest,
  kLinearWave,
  kGeostrophicAdjustment,
  kWilliamson2,
  kWilliamson6,
  kGalewsky
};
enum class DiffusionKind { kNone, kLaplacian, kBiharmonic };
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

struct ShallowWaterParameters {
  ShallowWaterTestCase test_case = ShallowWaterTestCase::kWilliamson2;
  ShallowWaterScheme scheme = ShallowWaterScheme::kRusanov;
  ReconstructionKind reconstruction = ReconstructionKind::kLinear;
  LimiterKind limiter = LimiterKind::kBarthJespersen;
  Real cfl = 0.5;
  Real mean_depth_m = 0.0;
  Real depth_floor_m = 0.0;
  DiffusionKind diffusion_kind = DiffusionKind::kNone;
  Real diffusion_coefficient = 0.0;
  Real flow_axis_x = 0.0;
  Real flow_axis_y = 0.0;
  Real flow_axis_z = 1.0;
  Real maximum_velocity_m_s = 0.0;
};

struct VerticalColumnParameters {
  VerticalTestCase test_case = VerticalTestCase::kIsothermal;
  Index levels = 0;
  std::vector<Real> a_half_pa;
  std::vector<Real> b_half;
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

struct DiagnosticsParameters {
  std::uint64_t interval_steps = 1;
};

struct ExperimentConfig {
  ExperimentKind kind = ExperimentKind::kOde;
  PlanetParameters planet;
  RunParameters run;
  OdeParameters ode;
  GridParameters grid{};
  TransportParameters transport{};
  ShallowWaterParameters shallow_water{};
  VerticalColumnParameters vertical{};
  DiagnosticsParameters diagnostics{};
  std::string output_directory;

  void validate() const;
};

[[nodiscard]] std::string_view experiment_kind_name(ExperimentKind kind) noexcept;
[[nodiscard]] std::string_view transport_scheme_name(TransportScheme scheme) noexcept;
[[nodiscard]] std::string_view limiter_name(LimiterKind limiter) noexcept;
[[nodiscard]] std::string_view transport_test_case_name(
    TransportTestCase test_case) noexcept;
[[nodiscard]] std::string_view initial_condition_name(
    InitialConditionKind initial_condition) noexcept;
[[nodiscard]] std::string_view shallow_water_scheme_name(
    ShallowWaterScheme scheme) noexcept;
[[nodiscard]] std::string_view reconstruction_name(
    ReconstructionKind reconstruction) noexcept;
[[nodiscard]] std::string_view shallow_water_test_case_name(
    ShallowWaterTestCase test_case) noexcept;
[[nodiscard]] std::string_view diffusion_kind_name(DiffusionKind kind) noexcept;
[[nodiscard]] std::string_view vertical_test_case_name(
    VerticalTestCase test_case) noexcept;
[[nodiscard]] std::string_view vertical_transport_scheme_name(
    VerticalTransportScheme scheme) noexcept;
[[nodiscard]] std::string_view vertical_limiter_name(
    VerticalLimiterKind limiter) noexcept;

[[nodiscard]] ExperimentConfig parse_experiment_config(std::istream& input);
[[nodiscard]] ExperimentConfig load_experiment_config(
    const std::filesystem::path& path);
void write_experiment_config(std::ostream& output, const ExperimentConfig& config);

}  // namespace mps
