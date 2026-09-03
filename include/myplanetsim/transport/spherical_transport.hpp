#pragma once

#include <cstdint>
#include <optional>
#include <span>
#include <vector>

#include "myplanetsim/config/experiment_config.hpp"
#include "myplanetsim/grid/cubed_sphere_grid.hpp"

namespace mps {

struct TransportState {
  Real time_s;
  std::uint64_t step;
  std::vector<Real> tracer;
};

struct TransportDiagnostics {
  Real initial_mass;
  Real final_mass;
  Real relative_mass_drift;
  Real l1_error;
  Real l2_error;
  Real linf_error;
  Real minimum;
  Real maximum;
  Real seam_rms_error;
  Real corner_rms_error;
  Real filament_preservation;
  Real unmixing;
  Real overshooting;
  Real real_mixing;
  std::uint64_t limiter_activations;
};

struct TransportResult {
  TransportState state;
  TransportDiagnostics diagnostics;
  bool reached_end_time;
};

[[nodiscard]] Vec3 solid_body_velocity(Vec3 unit_position, Real radius_m,
                                       Vec3 rotation_axis, Real angular_speed_rad_s);
[[nodiscard]] std::vector<Real> prescribed_edge_fluxes(const CubedSphereGrid& grid,
                                                       TransportTestCase test_case,
                                                       Real time_s, Real period_s,
                                                       Vec3 rotation_axis,
                                                       Real angular_speed_rad_s);
[[nodiscard]] Real initial_tracer_value(InitialConditionKind kind, Vec3 unit_position);
[[nodiscard]] std::vector<Real> make_initial_tracer(const CubedSphereGrid& grid,
                                                    InitialConditionKind kind);
[[nodiscard]] std::vector<Real> exact_solid_body_tracer(const CubedSphereGrid& grid,
                                                        InitialConditionKind kind,
                                                        Vec3 rotation_axis,
                                                        Real angular_speed_rad_s,
                                                        Real elapsed_time_s);
[[nodiscard]] Real tracer_mass(const CubedSphereGrid& grid,
                               std::span<const Real> tracer);
[[nodiscard]] Real stable_transport_time_step(const CubedSphereGrid& grid,
                                              std::span<const Real> oriented_edge_flux,
                                              Real cfl, Real maximum_time_step_s);

class SphericalTransportRhs {
 public:
  SphericalTransportRhs(const CubedSphereGrid& grid,
                        const TransportParameters& parameters, Real period_s);
  void operator()(Real time_s, std::span<const Real> tracer, std::span<Real> tendency);
  [[nodiscard]] std::uint64_t limiter_activations() const noexcept {
    return limiter_activations_;
  }

 private:
  const CubedSphereGrid& grid_;
  TransportParameters parameters_;
  Real period_s_;
  Vec3 axis_;
  std::uint64_t limiter_activations_ = 0;
};

[[nodiscard]] TransportResult run_spherical_transport(
    const ExperimentConfig& config,
    std::optional<TransportState> initial_state = std::nullopt,
    std::optional<std::uint64_t> stop_after_step = std::nullopt);

}  // namespace mps
