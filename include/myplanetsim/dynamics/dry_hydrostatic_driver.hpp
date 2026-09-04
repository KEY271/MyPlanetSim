#pragma once
#include <functional>

#include "myplanetsim/config/experiment_config.hpp"
#include "myplanetsim/dynamics/dry_hydrostatic_coupling.hpp"
#include "myplanetsim/grid/cubed_sphere_grid.hpp"
#include "myplanetsim/dynamics/surface_orography.hpp"
namespace mps {
struct DryHydrostaticRhs {
  std::vector<Real> surface_pressure_pa_s;
  DryHydrostaticTransportTendency tendency;
  Real horizontal_stable_time_step_s;
  Real maximum_continuity_residual_pa_s;
};
using DryHydrostaticObserver =
    std::function<void(const DryHydrostaticState&, const DryHydrostaticDerived&)>;
using DryHydrostaticCancel = std::function<bool()>;
class DryHydrostaticDriver {
 public:
  explicit DryHydrostaticDriver(ExperimentConfig config);
  [[nodiscard]] DryHydrostaticState initial_state() const;
  [[nodiscard]] DryHydrostaticDerived diagnose(const DryHydrostaticState&) const;
  [[nodiscard]] DryHydrostaticRhs rhs(const DryHydrostaticState&) const;
  void advance(DryHydrostaticState&, Real end_time_s,
               const DryHydrostaticObserver& observer = {},
               const DryHydrostaticCancel& cancel = {}) const;
  [[nodiscard]] const CubedSphereGrid& grid() const noexcept { return grid_; }
  [[nodiscard]] const SurfaceOrography& orography() const noexcept {
    return orography_;
  }

 private:
  ExperimentConfig config_;
  CubedSphereGrid grid_;
  AtmosphericHybridCoordinate coordinate_;
  SurfaceOrography orography_;
};
}  // namespace mps
