#pragma once
#include <functional>
#include <limits>
#include <optional>

#include "myplanetsim/config/experiment_config.hpp"
#include "myplanetsim/dynamics/dry_hydrostatic_coupling.hpp"
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
  Real horizontal_stable_time_step_s;
  Real vertical_stable_time_step_s;
  Real surface_stable_time_step_s = std::numeric_limits<Real>::infinity();
  Real maximum_continuity_residual_pa_s;
  HeldSuarezDiagnostics physics_diagnostics{};
  std::vector<Real> surface_temperature_k_s;
  SurfaceEnergyDiagnostics surface_diagnostics{};
  // Explicit horizontal diffusion, kept out of the advection, pressure-gradient,
  // Coriolis, and vertical-transport terms so its dissipation is attributable.
  Real diffusion_kinetic_energy_rate_w = 0.0;
};

struct DryHydrostaticStepDiagnostics {
  HeldSuarezDiagnostics physics_rates{};
  Real thermal_energy_contribution_j = 0.0;
  Real rayleigh_drag_energy_contribution_j = 0.0;
  SurfaceEnergyDiagnostics surface_rates{};
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
  [[nodiscard]] DryHydrostaticRhs rhs(const DryHydrostaticState&) const;
  void rhs(const DryHydrostaticState&, DryHydrostaticRhs& result) const;
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

 private:
  ExperimentConfig config_;
  CubedSphereGrid grid_;
  AtmosphericHybridCoordinate coordinate_;
  SurfaceOrography orography_;
  std::optional<SurfaceBoundary> surface_boundary_;
  mutable DryHydrostaticWorkspace workspace_;
};
}  // namespace mps
