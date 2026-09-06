#pragma once
#include <functional>
#include <limits>
#include <optional>

#include "myplanetsim/config/experiment_config.hpp"
#include "myplanetsim/dynamics/dry_hydrostatic_coupling.hpp"
#include "myplanetsim/dynamics/dry_hydrostatic_fast_modes.hpp"
#include "myplanetsim/dynamics/dry_hydrostatic_workspace.hpp"
#include "myplanetsim/dynamics/surface_boundary.hpp"
#include "myplanetsim/dynamics/surface_orography.hpp"
#include "myplanetsim/grid/cubed_sphere_grid.hpp"
#include "myplanetsim/physics/held_suarez.hpp"
#include "myplanetsim/physics/surface_energy_balance.hpp"
namespace mps {
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
  Real maximum_continuity_residual_pa_s;
  HeldSuarezDiagnostics physics_diagnostics{};
  std::vector<Real> surface_temperature_k_s;
  SurfaceEnergyDiagnostics surface_diagnostics{};
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

struct DryHydrostaticStepDiagnostics {
  HeldSuarezDiagnostics physics_rates{};
  Real thermal_energy_contribution_j = 0.0;
  Real rayleigh_drag_energy_contribution_j = 0.0;
  SurfaceEnergyDiagnostics surface_rates{};
  SurfaceEnergyBudget surface_budget{};
  Real diffusion_energy_contribution_j = 0.0;
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

 private:
  ExperimentConfig config_;
  CubedSphereGrid grid_;
  AtmosphericHybridCoordinate coordinate_;
  SurfaceOrography orography_;
  std::optional<SurfaceBoundary> surface_boundary_;
  std::optional<DryHydrostaticPressureReference> pressure_reference_;
  std::optional<DryHydrostaticReferenceColumn> semi_implicit_reference_column_;
  std::optional<DryHydrostaticVerticalModes> semi_implicit_vertical_modes_;
  mutable DryHydrostaticWorkspace workspace_;

  void rhs_with_components(const DryHydrostaticState&, DryHydrostaticRhs&,
                           DryHydrostaticRhsComponents*) const;
};
}  // namespace mps
