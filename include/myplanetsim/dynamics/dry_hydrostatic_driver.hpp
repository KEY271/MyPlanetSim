#pragma once
#include <array>
#include <functional>
#include <limits>
#include <optional>

#include "myplanetsim/config/experiment_config.hpp"
#include "myplanetsim/dynamics/dry_hydrostatic_coupling.hpp"
#include "myplanetsim/dynamics/dry_hydrostatic_fast_operator.hpp"
#include "myplanetsim/dynamics/dry_hydrostatic_semi_implicit.hpp"
#include "myplanetsim/dynamics/dry_hydrostatic_workspace.hpp"
#include "myplanetsim/dynamics/surface_boundary.hpp"
#include "myplanetsim/dynamics/surface_orography.hpp"
#include "myplanetsim/grid/cubed_sphere_grid.hpp"
#include "myplanetsim/physics/dry_mixing_coupling.hpp"
#include "myplanetsim/physics/gray_radiation_coupling.hpp"
#include "myplanetsim/physics/held_suarez.hpp"
#include "myplanetsim/physics/moist_physics_coupling.hpp"
#include "myplanetsim/physics/surface_energy_balance.hpp"
namespace mps {
// Exclusive RHS regions; reconstruction includes level packing and limiting.
struct DryHydrostaticRhsProfile {
  std::array<Real, 8> seconds{};
  std::size_t calls = 0;
};
struct DryHydrostaticRhs {
  std::vector<Real> surface_pressure_pa_s;
  DryHydrostaticTransportTendency tendency;
  // Legacy explicit limit. This remains min(fast-wave, diffusion) so the SSP-RK3
  // path is unchanged while Phase 10 can inspect the split constraints below.
  Real horizontal_stable_time_step_s = std::numeric_limits<Real>::infinity();
  Real horizontal_fast_wave_stable_time_step_s = std::numeric_limits<Real>::infinity();
  Real horizontal_advective_stable_time_step_s = std::numeric_limits<Real>::infinity();
  Real diffusion_stable_time_step_s = std::numeric_limits<Real>::infinity();
  Real vertical_stable_time_step_s = std::numeric_limits<Real>::infinity();
  Real surface_stable_time_step_s = std::numeric_limits<Real>::infinity();
  Real radiation_stable_time_step_s = std::numeric_limits<Real>::infinity();
  Real maximum_continuity_residual_pa_s;
  HeldSuarezDiagnostics physics_diagnostics{};
  std::vector<Real> surface_temperature_k_s;
  SurfaceEnergyDiagnostics surface_diagnostics{};
  GrayRadiationDiagnostics radiation_diagnostics{};
  std::size_t radiation_column_call_count = 0;
  Real radiation_wall_seconds = 0.0;
  // Explicit horizontal diffusion, kept out of the advection, pressure-gradient,
  // Coriolis, and vertical-transport terms so its dissipation is attributable.
  Real diffusion_kinetic_energy_rate_w = 0.0;
};

struct DryHydrostaticRhsTerm {
  std::vector<Real> surface_pressure_pa_s;
  DryHydrostaticTransportTendency tendency;
  std::vector<Real> surface_temperature_k_s;
};

struct DryHydrostaticRhsComponents {
  DryHydrostaticRhsTerm horizontal_transport;
  DryHydrostaticRhsTerm vertical_transport;
  DryHydrostaticRhsTerm pressure_gradient;
  DryHydrostaticRhsTerm coriolis;
  DryHydrostaticRhsTerm diffusion;
  DryHydrostaticRhsTerm physics;
};

struct DryConvectionDiagnostics {
  Real minimum_theta_difference_before_k = 0.0;
  Real minimum_theta_difference_after_k = 0.0;
  Real unstable_interface_fraction_before = 0.0;
  Real unstable_interface_fraction_after = 0.0;
  std::size_t adjusted_column_count = 0;
  std::size_t adjusted_layer_count = 0;
  std::size_t adjusted_block_count = 0;
  Real maximum_temperature_increment_k = 0.0;
  Real enthalpy_change_j = 0.0;
  Real dry_energy_attributed_change_j = 0.0;
};

// Categories are disjoint; discarded is a subset of iteration + final, not
// another category to add to total. Counts include failed attempts.
struct FullRhsEvaluationCounts {
  std::size_t initial = 0;
  std::size_t iteration = 0;
  std::size_t final = 0;
  std::size_t discarded = 0;
  [[nodiscard]] std::size_t total() const noexcept {
    return initial + iteration + final;
  }
};

struct DryHydrostaticStepDiagnostics {
  FullRhsEvaluationCounts full_rhs{};
  HeldSuarezDiagnostics physics_rates{};
  Real thermal_energy_contribution_j = 0.0;
  Real rayleigh_drag_energy_contribution_j = 0.0;
  SurfaceEnergyDiagnostics surface_rates{};
  SurfaceEnergyBudget surface_budget{};
  GrayRadiationDiagnostics radiation_rates{};
  GrayRadiationBudget radiation_budget{};
  DryMixingStepDiagnostics boundary_layer{};
  DryConvectionDiagnostics convection{};
  MoistPhysicsStepDiagnostics moisture{};
  Real diffusion_energy_contribution_j = 0.0;
  Real requested_time_step_s = 0.0;
  Real accepted_time_step_s = 0.0;
  Real radiation_stable_time_step_s = std::numeric_limits<Real>::infinity();
  Real advective_cfl = 0.0;
  Real implicit_wave_courant = 0.0;
  Real vertical_cfl = 0.0;
  std::size_t selected_implicit_modes = 0;
  std::size_t linear_iterations_total = 0;
  std::size_t linear_iterations_maximum = 0;
  Real linear_relative_residual_maximum = 0.0;
  std::size_t nonlinear_iterations = 0;
  Real nonlinear_relative_residual = 0.0;
  std::size_t retry_count = 0;
  std::size_t cfl_retry_count = 0;
  std::size_t invariant_retry_count = 0;
  std::size_t solver_retry_count = 0;
  std::size_t radiation_column_call_count = 0;
  std::size_t boundary_layer_column_call_count = 0;
  std::size_t convection_column_call_count = 0;
  std::size_t physics_substep_count = 0;
  std::size_t physics_retry_count = 0;
  Real radiation_wall_seconds = 0.0;
  Real wall_seconds_rhs = 0.0;
  Real wall_seconds_linear_solve = 0.0;
  Real wall_seconds_total = 0.0;
};
// The derived pointer is non-null at the configured diagnostic interval, initially,
// and at the final time. Step-integrated budgets are delivered on every accepted step.
using DryHydrostaticObserver =
    std::function<void(const DryHydrostaticState&, const DryHydrostaticDerived*,
                       const DryHydrostaticStepDiagnostics&)>;
using DryHydrostaticCancel = std::function<bool()>;
class DryHydrostaticDriver {
 public:
  explicit DryHydrostaticDriver(ExperimentConfig config);
  // Benchmark-only instrumentation, disabled by default; never serialized.
  void enable_rhs_profiling(bool enabled) const { profile_enabled_ = enabled; }
  void reset_rhs_profile() const { rhs_profile_ = {}; }
  [[nodiscard]] const DryHydrostaticRhsProfile& rhs_profile() const { return rhs_profile_; }
  [[nodiscard]] DryHydrostaticState initial_state() const;
  [[nodiscard]] DryHydrostaticDerived diagnose(const DryHydrostaticState&) const;
  [[nodiscard]] DryHydrostaticSources diagnose_sources(
      const DryHydrostaticDerived&) const;
  [[nodiscard]] DryHydrostaticRhs rhs(const DryHydrostaticState&) const;
  void rhs(const DryHydrostaticState&, DryHydrostaticRhs& result) const;
  [[nodiscard]] DryHydrostaticRhsComponents rhs_components(
      const DryHydrostaticState&) const;
  void advance(DryHydrostaticState&, Real end_time_s,
               const DryHydrostaticObserver& observer = {},
               const DryHydrostaticCancel& cancel = {}) const;
  [[nodiscard]] const CubedSphereGrid& grid() const noexcept { return grid_; }
  [[nodiscard]] const SurfaceOrography& orography() const noexcept {
    return orography_;
  }
  [[nodiscard]] const std::optional<SurfaceBoundary>& surface_boundary()
      const noexcept {
    return surface_boundary_;
  }
  [[nodiscard]] const std::optional<DryHydrostaticReferenceColumn>&
  semi_implicit_reference_column() const noexcept {
    return semi_implicit_reference_column_;
  }
  [[nodiscard]] const std::optional<DryHydrostaticVerticalModes>&
  semi_implicit_vertical_modes() const noexcept {
    return semi_implicit_vertical_modes_;
  }
  [[nodiscard]] const std::optional<DryHydrostaticFastOperator>&
  semi_implicit_fast_operator() const noexcept {
    return semi_implicit_fast_operator_;
  }

 private:
  mutable bool profile_enabled_ = false;
  mutable DryHydrostaticRhsProfile rhs_profile_{};
  ExperimentConfig config_;
  CubedSphereGrid grid_;
  AtmosphericHybridCoordinate coordinate_;
  SurfaceOrography orography_;
  std::optional<SurfaceBoundary> surface_boundary_;
  std::optional<DryHydrostaticPressureReference> pressure_reference_;
  // Mutable because the per-step reference update rebuilds them from the current state
  // inside the const `advance` path, like the reusable workspaces below.
  mutable std::optional<DryHydrostaticReferenceColumn> semi_implicit_reference_column_;
  mutable std::optional<DryHydrostaticVerticalModes> semi_implicit_vertical_modes_;
  mutable std::optional<DryHydrostaticFastOperator> semi_implicit_fast_operator_;
  mutable std::optional<DryHydrostaticExternalModeOperator>
      semi_implicit_external_operator_;
  mutable DryHydrostaticWorkspace workspace_;
  mutable DryHydrostaticSemiImplicitWorkspace semi_implicit_workspace_;

  void rhs_with_components(const DryHydrostaticState&, DryHydrostaticRhs&,
                           DryHydrostaticRhsComponents*,
                           bool compute_fast_wave_cfl = true) const;

  // Rebuilds the horizontally uniform reference column, fast operator, vertical modes
  // and external operator from the mass-weighted horizontal mean of `state`.
  void update_semi_implicit_reference(const DryHydrostaticState& state) const;
};
}  // namespace mps
