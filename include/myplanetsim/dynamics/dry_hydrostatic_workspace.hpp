#pragma once

#include <exception>

#include "myplanetsim/dynamics/dry_hydrostatic_coupling.hpp"
#include "myplanetsim/dynamics/dry_hydrostatic_reconstruction.hpp"
#include "myplanetsim/dynamics/dry_hydrostatic_sources.hpp"
#include "myplanetsim/physics/dry_convective_adjustment.hpp"
#include "myplanetsim/physics/dry_mixing_coupling.hpp"
#include "myplanetsim/physics/gray_radiation_coupling.hpp"
#include "myplanetsim/physics/held_suarez.hpp"
#include "myplanetsim/physics/moist_physics_coupling.hpp"
#include "myplanetsim/physics/surface_energy_balance.hpp"

namespace mps {

struct DryHydrostaticWorkspace {
  DryHydrostaticDerived derived;
  DryHydrostaticDiagnosisWorkspace diagnosis_workspace;

  DryHydrostaticTransportTendency horizontal_tendency;
  std::vector<DryHydrostaticEdgeFlux> edge_flux;
  std::vector<Real> edge_tracer_flux;
  std::vector<std::exception_ptr> edge_flux_failures;
  std::vector<Real> face_fast_wave_speed_length;
  std::vector<Real> face_advective_speed_length;
  DryHydrostaticReconstruction reconstruction;
  DryHydrostaticReconstructionWorkspace reconstruction_workspace;
  DryHydrostaticCoupling coupling;
  DryHydrostaticCouplingWorkspace coupling_workspace;
  DryHydrostaticSources sources;
  DryHydrostaticSourcesWorkspace sources_workspace;
  HeldSuarezTendency atmospheric_physics;
  HeldSuarezWorkspace atmospheric_physics_workspace;
  SurfaceEnergyTendency surface_physics;
  GrayRadiationTendency gray_radiation_physics;
  GrayRadiationCouplingWorkspace gray_radiation_workspace;
  DryConvectiveAdjustmentResult convective_adjustment;
  DryConvectiveAdjustmentWorkspace convective_adjustment_workspace;
  DryMixingStepDiagnostics dry_mixing_diagnostics;
  DryMixingCouplingWorkspace dry_mixing_workspace;
  MoistPhysicsStepDiagnostics moist_diagnostics;
  MoistPhysicsCouplingWorkspace moist_workspace;
};

}  // namespace mps
