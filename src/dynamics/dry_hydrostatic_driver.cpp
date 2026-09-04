#include "myplanetsim/dynamics/dry_hydrostatic_driver.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "myplanetsim/dynamics/dry_hydrostatic_benchmarks.hpp"
#include "myplanetsim/dynamics/dry_hydrostatic_diffusion.hpp"
#include "myplanetsim/dynamics/dry_hydrostatic_flux.hpp"
#include "myplanetsim/dynamics/dry_hydrostatic_reconstruction.hpp"
#include "myplanetsim/dynamics/dry_hydrostatic_sources.hpp"
#include "myplanetsim/dynamics/shallow_water_diffusion.hpp"
#include "myplanetsim/physics/held_suarez.hpp"
#include "myplanetsim/physics/planetary_newtonian.hpp"
#include "myplanetsim/physics/surface_energy_balance.hpp"
namespace mps {
namespace {

[[nodiscard]] DryHydrostaticState euler_update(const DryHydrostaticState& base,
                                               const DryHydrostaticRhs& tendency,
                                               const Real scale) {
  DryHydrostaticState result = base;
  for (std::size_t cell = 0; cell < result.surface_pressure_pa.size(); ++cell) {
    result.surface_pressure_pa[cell] += scale * tendency.surface_pressure_pa_s[cell];
  }
  for (std::size_t n = 0; n < result.horizontal_momentum_mass_kg_m_s.size(); ++n) {
    result.horizontal_momentum_mass_kg_m_s[n] =
        result.horizontal_momentum_mass_kg_m_s[n] +
        scale * tendency.tendency.momentum[n];
    result.potential_temperature_mass_k_kg_m2[n] +=
        scale * tendency.tendency.potential_temperature_mass[n];
    result.tracer_mass_kg_m2[n] += scale * tendency.tendency.tracer_mass[n];
  }
  for (std::size_t cell = 0; cell < result.surface_temperature_k.size(); ++cell)
    result.surface_temperature_k[cell] +=
        scale * tendency.surface_temperature_k_s[cell];
  return result;
}

[[nodiscard]] DryHydrostaticState convex_update(const DryHydrostaticState& initial,
                                                const Real initial_weight,
                                                const DryHydrostaticState& stage,
                                                const DryHydrostaticRhs& tendency,
                                                const Real stage_weight,
                                                const Real time_step_s) {
  DryHydrostaticState result = initial;
  for (std::size_t cell = 0; cell < result.surface_pressure_pa.size(); ++cell) {
    result.surface_pressure_pa[cell] =
        initial_weight * initial.surface_pressure_pa[cell] +
        stage_weight * (stage.surface_pressure_pa[cell] +
                        time_step_s * tendency.surface_pressure_pa_s[cell]);
  }
  for (std::size_t n = 0; n < result.horizontal_momentum_mass_kg_m_s.size(); ++n) {
    result.horizontal_momentum_mass_kg_m_s[n] =
        initial_weight * initial.horizontal_momentum_mass_kg_m_s[n] +
        stage_weight * (stage.horizontal_momentum_mass_kg_m_s[n] +
                        time_step_s * tendency.tendency.momentum[n]);
    result.potential_temperature_mass_k_kg_m2[n] =
        initial_weight * initial.potential_temperature_mass_k_kg_m2[n] +
        stage_weight * (stage.potential_temperature_mass_k_kg_m2[n] +
                        time_step_s * tendency.tendency.potential_temperature_mass[n]);
    result.tracer_mass_kg_m2[n] =
        initial_weight * initial.tracer_mass_kg_m2[n] +
        stage_weight * (stage.tracer_mass_kg_m2[n] +
                        time_step_s * tendency.tendency.tracer_mass[n]);
  }
  for (std::size_t cell = 0; cell < result.surface_temperature_k.size(); ++cell)
    result.surface_temperature_k[cell] =
        initial_weight * initial.surface_temperature_k[cell] +
        stage_weight * (stage.surface_temperature_k[cell] +
                        time_step_s * tendency.surface_temperature_k_s[cell]);
  return result;
}

[[nodiscard]] bool pressure_is_in_range(const ExperimentConfig& config,
                                        const DryHydrostaticState& state) {
  return std::ranges::all_of(state.surface_pressure_pa, [&](const Real pressure) {
    return std::isfinite(pressure) &&
           pressure >= config.vertical.minimum_surface_pressure_pa &&
           pressure <= config.vertical.maximum_surface_pressure_pa;
  });
}

[[nodiscard]] Real pressure_limited_time_step(const ExperimentConfig& config,
                                              const DryHydrostaticState& state,
                                              const DryHydrostaticRhs& tendency,
                                              const Real maximum_time_step_s) {
  Real result = maximum_time_step_s;
  for (std::size_t cell = 0; cell < state.surface_pressure_pa.size(); ++cell) {
    const Real rate = tendency.surface_pressure_pa_s[cell];
    if (rate == 0.0) {
      continue;
    }
    const Real available = rate > 0.0 ? config.vertical.maximum_surface_pressure_pa -
                                            state.surface_pressure_pa[cell]
                                      : state.surface_pressure_pa[cell] -
                                            config.vertical.minimum_surface_pressure_pa;
    if (!(available > 0.0)) {
      throw std::runtime_error(
          "dry dynamics drives surface pressure outside the configured range");
    }
    const Real limit = available / std::abs(rate);
    if (limit < result) {
      result = std::nextafter(limit, 0.0);
    }
  }
  return result;
}

void halve_time_step(const DryHydrostaticState& state, Real& time_step_s) {
  time_step_s *= 0.5;
  if (!(time_step_s > 0.0) || state.time_s + time_step_s == state.time_s) {
    throw std::runtime_error(
        "dry stage pressure or CFL constraint produced a zero time step");
  }
}

[[nodiscard]] Real stable_time_step(const DryHydrostaticRhs& rhs) {
  return std::min(rhs.horizontal_stable_time_step_s, rhs.vertical_stable_time_step_s);
}

[[nodiscard]] bool surface_temperature_is_positive(const DryHydrostaticState& state) {
  return std::ranges::all_of(state.surface_temperature_k, [](const Real value) {
    return value > 0.0 && std::isfinite(value);
  });
}

}  // namespace

DryHydrostaticDriver::DryHydrostaticDriver(ExperimentConfig c)
    : config_(std::move(c)),
      grid_(config_.grid.cells_per_panel, config_.planet.radius_m),
      coordinate_({config_.vertical.a_half_pa, config_.vertical.b_half},
                  config_.vertical.minimum_surface_pressure_pa,
                  config_.vertical.maximum_surface_pressure_pa,
                  config_.vertical.minimum_pressure_thickness_pa),
      orography_(make_surface_orography(config_.orography, grid_,
                                        config_.planet.gravity_m_s2,
                                        config_.source_directory)) {
  if (config_.kind != ExperimentKind::kDryHydrostatic)
    throw std::invalid_argument("dry driver requires dry_hydrostatic config");
  if (config_.surface.has_value()) {
    surface_boundary_ =
        make_surface_boundary(*config_.surface, config_.orography, grid_,
                              config_.planet.gravity_m_s2, config_.source_directory);
    orography_ = SurfaceOrography(
        std::vector<Real>(surface_boundary_->surface_geopotential_m2_s2().begin(),
                          surface_boundary_->surface_geopotential_m2_s2().end()),
        surface_boundary_->source_fingerprint());
  }
}
DryHydrostaticState DryHydrostaticDriver::initial_state() const {
  auto state =
      initialize_dry_hydrostatic_benchmark(config_, grid_, coordinate_, orography_);
  if (config_.physics.kind == PhysicsKind::kSurfaceEnergyBalance)
    state.surface_temperature_k.assign(grid_.cell_count(),
                                       config_.surface->initial_temperature_k);
  return state;
}
DryHydrostaticDerived DryHydrostaticDriver::diagnose(
    const DryHydrostaticState& s) const {
  return diagnose_dry_hydrostatic_state(s, coordinate_, config_.planet,
                                        orography_.surface_geopotential_m2_s2());
}
DryHydrostaticRhs DryHydrostaticDriver::rhs(const DryHydrostaticState& s) const {
  auto d = diagnose(s);
  auto C = d.cells, K = d.levels;
  DryHydrostaticTransportTendency h;
  h.air_mass.assign(C * K, 0);
  h.momentum.assign(C * K, {});
  h.potential_temperature_mass.assign(C * K, 0);
  h.tracer_mass.assign(C * K, 0);
  Real dt = std::numeric_limits<Real>::infinity();
  const auto reconstructed = reconstruct_dry_hydrostatic_face_states(
      grid_, d, config_.dry_hydrostatic.reconstruction,
      config_.dry_hydrostatic.limiter);
  for (const auto& e : grid_.edges()) {
    auto l = grid_.cell_index(e.left_cell);
    auto r = grid_.cell_index(e.right_cell);
    auto basis = edge_tangent_basis(e);
    for (std::size_t k = 0; k < K; ++k) {
      const auto& face = reconstructed.at(e.id, k);
      auto f = rusanov_dry_hydrostatic_flux(face.left, face.right, basis,
                                            config_.planet.gas_constant_j_kg_k,
                                            config_.planet.heat_capacity_cp_j_kg_k);
      auto add = [&](std::size_t c, Real sign) {
        auto n = dry_hydrostatic_offset(c, k, K);
        auto scale = sign * e.length_m / grid_.cells()[c].area_m2;
        h.air_mass[n] += scale * f.air_mass_kg_m_s;
        h.momentum[n] = h.momentum[n] + scale * f.momentum_kg_s2;
        h.potential_temperature_mass[n] +=
            scale * f.potential_temperature_mass_k_kg_m_s;
        h.tracer_mass[n] += scale * f.tracer_mass_kg_m_s;
      };
      add(l, -1);
      add(r, 1);
      dt = std::min(dt,
                    config_.dry_hydrostatic.cfl *
                        std::min(grid_.cells()[l].area_m2, grid_.cells()[r].area_m2) /
                        (e.length_m * f.maximum_wave_speed_m_s));
    }
  }
  auto coupled = couple_dry_hydrostatic_columns(
      s, d, h, coordinate_.coefficients().b_half, config_.planet.gravity_m_s2,
      config_.vertical.transport_scheme, config_.vertical.limiter);
  Real vertical_dt = config_.run.time_step_s;
  for (std::size_t c = 0; c < C; ++c) {
    const auto begin = c * K;
    const auto interface_begin = c * (K + 1);
    const std::vector<Real> air_mass(
        d.air_mass_kg_m2.begin() + static_cast<std::ptrdiff_t>(begin),
        d.air_mass_kg_m2.begin() + static_cast<std::ptrdiff_t>(begin + K));
    const std::vector<Real> interface_flux(
        coupled.interface_mass_flux_kg_m2_s.begin() +
            static_cast<std::ptrdiff_t>(interface_begin),
        coupled.interface_mass_flux_kg_m2_s.begin() +
            static_cast<std::ptrdiff_t>(interface_begin + K + 1));
    const std::vector<Real> horizontal_air_mass(
        h.air_mass.begin() + static_cast<std::ptrdiff_t>(begin),
        h.air_mass.begin() + static_cast<std::ptrdiff_t>(begin + K));
    vertical_dt = std::min(
        vertical_dt,
        vertical_stable_time_step(air_mass, interface_flux, horizontal_air_mass,
                                  config_.vertical.cfl, config_.run.time_step_s));
  }
  auto sources = dry_hydrostatic_sources(grid_, d, config_.planet);
  for (std::size_t n = 0; n < C * K; ++n)
    coupled.tendency.momentum[n] = coupled.tendency.momentum[n] +
                                   sources.pressure_gradient_kg_m_s2[n] +
                                   sources.coriolis_kg_m_s2[n];
  Real diffusion_rate = 0.0;
  if (config_.dry_hydrostatic.diffusion_kind != DiffusionKind::kNone) {
    const auto diffusion = dry_hydrostatic_diffusion_tendency(
        grid_, d, config_.dry_hydrostatic.diffusion_kind,
        config_.dry_hydrostatic.diffusion_coefficient);
    for (std::size_t n = 0; n < C * K; ++n) {
      coupled.tendency.momentum[n] =
          coupled.tendency.momentum[n] + diffusion.momentum[n];
      coupled.tendency.potential_temperature_mass[n] +=
          diffusion.potential_temperature_mass[n];
      coupled.tendency.tracer_mass[n] += diffusion.tracer_mass[n];
    }
    diffusion_rate = diffusion.kinetic_energy_rate_w;
    dt = std::min(
        dt, stable_diffusion_time_step(grid_, config_.dry_hydrostatic.diffusion_kind,
                                       config_.dry_hydrostatic.diffusion_coefficient,
                                       config_.run.time_step_s));
  }
  HeldSuarezDiagnostics physics_diagnostics{};
  SurfaceEnergyDiagnostics surface_diagnostics{};
  std::vector<Real> surface_temperature_rate;
  if (config_.physics.kind == PhysicsKind::kHeldSuarez) {
    auto physics = held_suarez_tendency(grid_, coordinate_, d, s.surface_pressure_pa,
                                        config_.planet);
    for (std::size_t n = 0; n < C * K; ++n) {
      coupled.tendency.momentum[n] =
          coupled.tendency.momentum[n] + physics.horizontal_momentum_mass_kg_m_s2[n];
      coupled.tendency.potential_temperature_mass[n] +=
          physics.potential_temperature_mass_k_kg_m2_s[n];
    }
    physics_diagnostics = physics.diagnostics;
  } else if (config_.physics.kind == PhysicsKind::kSurfaceEnergyBalance) {
    const auto orbit_state =
        evaluate_orbit(*config_.orbit, config_.planet.rotation_rate_rad_s,
                       s.time_s - config_.run.start_time_s);
    auto physics = surface_energy_tendency(
        grid_, *surface_boundary_, s.surface_temperature_k, d, s.surface_pressure_pa,
        config_.planet, *config_.surface, orbit_state);
    for (std::size_t n = 0; n < C * K; ++n) {
      coupled.tendency.momentum[n] =
          coupled.tendency.momentum[n] + physics.horizontal_momentum_mass_kg_m_s2[n];
      coupled.tendency.potential_temperature_mass[n] +=
          physics.potential_temperature_mass_k_kg_m2_s[n];
    }
    surface_temperature_rate = std::move(physics.surface_temperature_k_s);
    surface_diagnostics = physics.diagnostics;
  } else if (config_.physics.kind == PhysicsKind::kPlanetaryNewtonian) {
    std::optional<OrbitState> orbit_state;
    if (config_.physics.geometry == ForcingGeometry::kSubstellar)
      orbit_state = evaluate_orbit(*config_.orbit, config_.planet.rotation_rate_rad_s,
                                   s.time_s - config_.run.start_time_s);
    auto physics = planetary_newtonian_tendency(
        grid_, coordinate_, d, s.surface_pressure_pa, config_.planet,
        config_.physics.geometry, orbit_state ? &*orbit_state : nullptr);
    for (std::size_t n = 0; n < C * K; ++n) {
      coupled.tendency.momentum[n] =
          coupled.tendency.momentum[n] + physics.horizontal_momentum_mass_kg_m_s2[n];
      coupled.tendency.potential_temperature_mass[n] +=
          physics.potential_temperature_mass_k_kg_m2_s[n];
    }
    physics_diagnostics = physics.diagnostics;
  }
  return {.surface_pressure_pa_s = std::move(coupled.surface_pressure_pa_s),
          .tendency = std::move(coupled.tendency),
          .horizontal_stable_time_step_s = dt,
          .vertical_stable_time_step_s = vertical_dt,
          .maximum_continuity_residual_pa_s = coupled.maximum_continuity_residual_pa_s,
          .physics_diagnostics = physics_diagnostics,
          .surface_temperature_k_s = std::move(surface_temperature_rate),
          .surface_diagnostics = surface_diagnostics,
          .diffusion_kinetic_energy_rate_w = diffusion_rate};
}
void DryHydrostaticDriver::advance(DryHydrostaticState& s, const Real end,
                                   const DryHydrostaticObserver& obs,
                                   const DryHydrostaticCancel& cancel) const {
  const auto observe = [&](const DryHydrostaticState& state,
                           const DryHydrostaticStepDiagnostics& step) {
    if (!obs) return;
    const auto derived = diagnose(state);
    auto sampled = step;
    if (config_.physics.kind == PhysicsKind::kHeldSuarez) {
      sampled.physics_rates =
          held_suarez_tendency(grid_, coordinate_, derived, state.surface_pressure_pa,
                               config_.planet)
              .diagnostics;
    } else if (config_.physics.kind == PhysicsKind::kPlanetaryNewtonian) {
      std::optional<OrbitState> orbit_state;
      if (config_.physics.geometry == ForcingGeometry::kSubstellar)
        orbit_state = evaluate_orbit(*config_.orbit, config_.planet.rotation_rate_rad_s,
                                     state.time_s - config_.run.start_time_s);
      sampled.physics_rates =
          planetary_newtonian_tendency(
              grid_, coordinate_, derived, state.surface_pressure_pa, config_.planet,
              config_.physics.geometry, orbit_state ? &*orbit_state : nullptr)
              .diagnostics;
    } else if (config_.physics.kind == PhysicsKind::kSurfaceEnergyBalance) {
      const auto orbit_state =
          evaluate_orbit(*config_.orbit, config_.planet.rotation_rate_rad_s,
                         state.time_s - config_.run.start_time_s);
      sampled.surface_rates =
          surface_energy_tendency(
              grid_, *surface_boundary_, state.surface_temperature_k, derived,
              state.surface_pressure_pa, config_.planet, *config_.surface, orbit_state)
              .diagnostics;
    }
    obs(state, derived, sampled);
  };
  observe(s, {});
  std::vector<Vec3> centres;
  centres.reserve(grid_.cell_count());
  for (const auto& cell : grid_.cells()) centres.push_back(cell.center);
  const auto project_and_validate = [&](DryHydrostaticState& stage) {
    for (std::size_t n = 0; n < stage.horizontal_momentum_mass_kg_m_s.size(); ++n) {
      const auto cell = n / coordinate_.levels();
      stage.horizontal_momentum_mass_kg_m_s[n] = project_tangent(
          stage.horizontal_momentum_mass_kg_m_s[n], grid_.cells()[cell].center);
    }
    validate_dry_hydrostatic_state(stage, diagnose(stage), centres,
                                   config_.vertical.minimum_surface_pressure_pa,
                                   config_.vertical.maximum_surface_pressure_pa,
                                   config_.vertical.temperature_floor_k);
  };
  while (s.time_s < end) {
    if (cancel && cancel()) return;
    const DryHydrostaticState initial = s;
    const auto rhs1 = rhs(initial);
    Real dt = std::min(
        {config_.run.time_step_s, stable_time_step(rhs1), end - initial.time_s});
    dt = pressure_limited_time_step(config_, initial, rhs1, dt);
    if (!(dt > 0.0) || !std::isfinite(dt)) {
      throw std::runtime_error("dry driver stable time step is invalid");
    }

    while (true) {
      auto stage1 = euler_update(initial, rhs1, dt);
      stage1.time_s = initial.time_s + dt;
      if (!pressure_is_in_range(config_, stage1) ||
          !surface_temperature_is_positive(stage1)) {
        halve_time_step(initial, dt);
        continue;
      }
      project_and_validate(stage1);
      const auto rhs2 = rhs(stage1);
      if (dt > stable_time_step(rhs2)) {
        halve_time_step(initial, dt);
        continue;
      }

      auto stage2 = convex_update(initial, 0.75, stage1, rhs2, 0.25, dt);
      stage2.time_s = initial.time_s + 0.5 * dt;
      if (!pressure_is_in_range(config_, stage2) ||
          !surface_temperature_is_positive(stage2)) {
        halve_time_step(initial, dt);
        continue;
      }
      project_and_validate(stage2);
      const auto rhs3 = rhs(stage2);
      if (dt > stable_time_step(rhs3)) {
        halve_time_step(initial, dt);
        continue;
      }

      auto next = convex_update(initial, 1.0 / 3.0, stage2, rhs3, 2.0 / 3.0, dt);
      next.time_s = initial.time_s + dt;
      next.step = initial.step + 1;
      if (!pressure_is_in_range(config_, next) ||
          !surface_temperature_is_positive(next)) {
        halve_time_step(initial, dt);
        continue;
      }
      project_and_validate(next);
      s = std::move(next);
      observe(s, {.thermal_energy_contribution_j =
                      dt * (rhs1.physics_diagnostics.thermal_energy_rate_w / 6.0 +
                            rhs2.physics_diagnostics.thermal_energy_rate_w / 6.0 +
                            2.0 * rhs3.physics_diagnostics.thermal_energy_rate_w / 3.0),
                  .rayleigh_drag_energy_contribution_j =
                      dt * (rhs1.physics_diagnostics.rayleigh_drag_work_w / 6.0 +
                            rhs2.physics_diagnostics.rayleigh_drag_work_w / 6.0 +
                            2.0 * rhs3.physics_diagnostics.rayleigh_drag_work_w / 3.0),
                  .diffusion_energy_contribution_j =
                      dt * (rhs1.diffusion_kinetic_energy_rate_w / 6.0 +
                            rhs2.diffusion_kinetic_energy_rate_w / 6.0 +
                            2.0 * rhs3.diffusion_kinetic_energy_rate_w / 3.0)});
      break;
    }
  }
}
}  // namespace mps
