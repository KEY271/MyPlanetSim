#include "myplanetsim/dynamics/dry_hydrostatic_driver.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>

#include "myplanetsim/dynamics/dry_hydrostatic_benchmarks.hpp"
#include "myplanetsim/dynamics/dry_hydrostatic_diffusion.hpp"
#include "myplanetsim/dynamics/dry_hydrostatic_flux.hpp"
#include "myplanetsim/dynamics/dry_hydrostatic_reconstruction.hpp"
#include "myplanetsim/dynamics/dry_hydrostatic_sources.hpp"
#include "myplanetsim/dynamics/shallow_water_diffusion.hpp"
#include "myplanetsim/physics/gray_radiation_coupling.hpp"
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
  return std::min({rhs.horizontal_stable_time_step_s, rhs.vertical_stable_time_step_s,
                   rhs.surface_stable_time_step_s, rhs.radiation_stable_time_step_s});
}

void resize_zero_rhs_term(DryHydrostaticRhsTerm& term, const std::size_t cells,
                          const std::size_t levels) {
  const auto cell_levels = cells * levels;
  term.surface_pressure_pa_s.assign(cells, 0.0);
  term.tendency.air_mass.assign(cell_levels, 0.0);
  term.tendency.momentum.assign(cell_levels, {});
  term.tendency.potential_temperature_mass.assign(cell_levels, 0.0);
  term.tendency.tracer_mass.assign(cell_levels, 0.0);
  term.surface_temperature_k_s.assign(cells, 0.0);
}

void resize_zero_rhs_components(DryHydrostaticRhsComponents& components,
                                const std::size_t cells, const std::size_t levels) {
  resize_zero_rhs_term(components.horizontal_transport, cells, levels);
  resize_zero_rhs_term(components.vertical_transport, cells, levels);
  resize_zero_rhs_term(components.pressure_gradient, cells, levels);
  resize_zero_rhs_term(components.coriolis, cells, levels);
  resize_zero_rhs_term(components.diffusion, cells, levels);
  resize_zero_rhs_term(components.physics, cells, levels);
}

[[nodiscard]] bool surface_temperature_is_positive(const DryHydrostaticState& state) {
  return std::ranges::all_of(state.surface_temperature_k, [](const Real value) {
    return value > 0.0 && std::isfinite(value);
  });
}

[[nodiscard]] DryHydrostaticState crank_nicolson_residual(
    const DryHydrostaticState& initial, const DryHydrostaticState& candidate,
    const DryHydrostaticRhs& initial_rhs, const DryHydrostaticRhs& candidate_rhs,
    const Real time_step_s, const Real implicit_weight,
    const std::span<const Vec3> cell_centres, const std::size_t levels) {
  DryHydrostaticState residual = candidate;
  const Real initial_weight = 1.0 - implicit_weight;
  for (std::size_t cell = 0; cell < residual.surface_pressure_pa.size(); ++cell) {
    residual.surface_pressure_pa[cell] =
        candidate.surface_pressure_pa[cell] - initial.surface_pressure_pa[cell] -
        time_step_s * (initial_weight * initial_rhs.surface_pressure_pa_s[cell] +
                       implicit_weight * candidate_rhs.surface_pressure_pa_s[cell]);
  }
  for (std::size_t n = 0; n < residual.horizontal_momentum_mass_kg_m_s.size(); ++n) {
    residual.horizontal_momentum_mass_kg_m_s[n] =
        candidate.horizontal_momentum_mass_kg_m_s[n] -
        initial.horizontal_momentum_mass_kg_m_s[n] -
        time_step_s * (initial_weight * initial_rhs.tendency.momentum[n] +
                       implicit_weight * candidate_rhs.tendency.momentum[n]);
    residual.horizontal_momentum_mass_kg_m_s[n] = project_tangent(
        residual.horizontal_momentum_mass_kg_m_s[n], cell_centres[n / levels]);
    residual.potential_temperature_mass_k_kg_m2[n] =
        candidate.potential_temperature_mass_k_kg_m2[n] -
        initial.potential_temperature_mass_k_kg_m2[n] -
        time_step_s *
            (initial_weight * initial_rhs.tendency.potential_temperature_mass[n] +
             implicit_weight * candidate_rhs.tendency.potential_temperature_mass[n]);
    residual.tracer_mass_kg_m2[n] =
        candidate.tracer_mass_kg_m2[n] - initial.tracer_mass_kg_m2[n] -
        time_step_s * (initial_weight * initial_rhs.tendency.tracer_mass[n] +
                       implicit_weight * candidate_rhs.tendency.tracer_mass[n]);
  }
  for (std::size_t cell = 0; cell < residual.surface_temperature_k.size(); ++cell) {
    residual.surface_temperature_k[cell] =
        candidate.surface_temperature_k[cell] - initial.surface_temperature_k[cell] -
        time_step_s * (initial_weight * initial_rhs.surface_temperature_k_s[cell] +
                       implicit_weight * candidate_rhs.surface_temperature_k_s[cell]);
  }
  return residual;
}

struct ScaledCrankNicolsonResidual {
  Real norm = 0.0;
  Real rms_norm = 0.0;
  const char* component = "none";
};

[[nodiscard]] ScaledCrankNicolsonResidual scaled_crank_nicolson_residual(
    const DryHydrostaticState& residual, const DryHydrostaticState& initial,
    const DryHydrostaticReferenceColumn& reference,
    const DryHydrostaticExternalModeOperator& external) {
  ScaledCrankNicolsonResidual result;
  Real sum_of_squares = 0.0;
  std::size_t sample_count = 0;
  const auto consider = [&](const Real value, const char* const component) {
    if (value > result.norm) {
      result.norm = value;
      result.component = component;
    }
    sum_of_squares += value * value;
    ++sample_count;
  };
  for (std::size_t cell = 0; cell < residual.surface_pressure_pa.size(); ++cell) {
    consider(std::abs(residual.surface_pressure_pa[cell]) /
                 std::max({1.0, reference.surface_pressure_pa,
                           std::abs(initial.surface_pressure_pa[cell])}),
             "surface_pressure");
  }
  const auto levels = reference.geometry.air_mass_kg_m2.size();
  for (std::size_t n = 0; n < residual.horizontal_momentum_mass_kg_m_s.size(); ++n) {
    const auto level = n % levels;
    const Real momentum_scale = std::max(
        1.0, reference.geometry.air_mass_kg_m2[level] * external.phase_speed_m_s);
    consider(
        norm(residual.horizontal_momentum_mass_kg_m_s[n]) /
            std::max(momentum_scale, norm(initial.horizontal_momentum_mass_kg_m_s[n])),
        "momentum");
    consider(std::abs(residual.potential_temperature_mass_k_kg_m2[n]) /
                 std::max({1.0, reference.potential_temperature_mass_k_kg_m2[level],
                           std::abs(initial.potential_temperature_mass_k_kg_m2[n])}),
             "potential_temperature_mass");
    consider(std::abs(residual.tracer_mass_kg_m2[n]) /
                 std::max(1.0, std::abs(initial.tracer_mass_kg_m2[n])),
             "tracer_mass");
  }
  for (std::size_t cell = 0; cell < residual.surface_temperature_k.size(); ++cell) {
    consider(std::abs(residual.surface_temperature_k[cell]) /
                 std::max({1.0, reference.temperature_k,
                           std::abs(initial.surface_temperature_k[cell])}),
             "surface_temperature");
  }
  result.rms_norm = std::sqrt(sum_of_squares / static_cast<Real>(sample_count));
  return result;
}

[[nodiscard]] DryHydrostaticFastPerturbation negative_fast_residual(
    const DryHydrostaticState& residual) {
  DryHydrostaticFastPerturbation result{
      .surface_pressure_pa = residual.surface_pressure_pa,
      .horizontal_momentum_mass_kg_m_s = residual.horizontal_momentum_mass_kg_m_s,
      .potential_temperature_mass_k_kg_m2 =
          residual.potential_temperature_mass_k_kg_m2};
  for (auto& value : result.surface_pressure_pa) value = -value;
  for (auto& value : result.horizontal_momentum_mass_kg_m_s) value = -value;
  for (auto& value : result.potential_temperature_mass_k_kg_m2) value = -value;
  return result;
}

void add_semi_implicit_correction(DryHydrostaticState& state,
                                  const DryHydrostaticFastPerturbation& fast_correction,
                                  const DryHydrostaticState& residual) {
  for (std::size_t cell = 0; cell < state.surface_pressure_pa.size(); ++cell)
    state.surface_pressure_pa[cell] += fast_correction.surface_pressure_pa[cell];
  for (std::size_t n = 0; n < state.horizontal_momentum_mass_kg_m_s.size(); ++n) {
    state.horizontal_momentum_mass_kg_m_s[n] =
        state.horizontal_momentum_mass_kg_m_s[n] +
        fast_correction.horizontal_momentum_mass_kg_m_s[n];
    state.potential_temperature_mass_k_kg_m2[n] +=
        fast_correction.potential_temperature_mass_k_kg_m2[n];
    // Tracer and surface reservoirs are outside L_ref, so their quasi-Newton
    // correction is the identity solve.
    state.tracer_mass_kg_m2[n] -= residual.tracer_mass_kg_m2[n];
  }
  for (std::size_t cell = 0; cell < state.surface_temperature_k.size(); ++cell)
    state.surface_temperature_k[cell] -= residual.surface_temperature_k[cell];
}

[[nodiscard]] Real semi_implicit_explicit_limit(const ExperimentConfig& config,
                                                const DryHydrostaticRhs& rhs) {
  const Real advective = rhs.horizontal_advective_stable_time_step_s *
                         config.dry_hydrostatic.advective_cfl /
                         config.dry_hydrostatic.cfl;
  return std::min({advective, rhs.diffusion_stable_time_step_s,
                   rhs.vertical_stable_time_step_s, rhs.surface_stable_time_step_s,
                   rhs.radiation_stable_time_step_s});
}

[[nodiscard]] Real courant_from_stable_time_step(const Real time_step_s,
                                                 const Real configured_cfl,
                                                 const Real stable_time_step_s) {
  if (std::isinf(stable_time_step_s)) return 0.0;
  if (!(stable_time_step_s > 0.0) || !std::isfinite(stable_time_step_s))
    throw std::runtime_error("stable time step is invalid for Courant diagnostics");
  return time_step_s * configured_cfl / stable_time_step_s;
}

[[nodiscard]] std::string nonlinear_failure_message(const Real initial_residual,
                                                    const Real residual,
                                                    const Real rms_residual,
                                                    const char* const component,
                                                    const Real time_step_s) {
  std::ostringstream message;
  message << std::setprecision(17)
          << "Crank-Nicolson iteration diverged: initial_residual=" << initial_residual
          << ", residual=" << residual << ", rms_residual=" << rms_residual
          << ", component=" << component << ", dt_s=" << time_step_s;
  return message.str();
}

[[nodiscard]] std::string modal_failure_message(const Real equation_residual) {
  std::ostringstream message;
  message << std::setprecision(17)
          << "modal linear solve did not converge: equation_residual="
          << equation_residual;
  return message.str();
}

[[nodiscard]] GrayRadiationDiagnostics weighted_radiation_diagnostics(
    const GrayRadiationDiagnostics& first, const Real first_weight,
    const GrayRadiationDiagnostics& second, const Real second_weight,
    const GrayRadiationDiagnostics& third = {}, const Real third_weight = 0.0) {
  const auto weighted = [&](const Real GrayRadiationDiagnostics::* member) {
    return first_weight * first.*member + second_weight * second.*member +
           third_weight * third.*member;
  };
  return {
      .toa_incoming_shortwave_power_w =
          weighted(&GrayRadiationDiagnostics::toa_incoming_shortwave_power_w),
      .toa_reflected_shortwave_power_w =
          weighted(&GrayRadiationDiagnostics::toa_reflected_shortwave_power_w),
      .toa_outgoing_longwave_power_w =
          weighted(&GrayRadiationDiagnostics::toa_outgoing_longwave_power_w),
      .toa_net_upward_power_w =
          weighted(&GrayRadiationDiagnostics::toa_net_upward_power_w),
      .surface_down_shortwave_power_w =
          weighted(&GrayRadiationDiagnostics::surface_down_shortwave_power_w),
      .surface_up_shortwave_power_w =
          weighted(&GrayRadiationDiagnostics::surface_up_shortwave_power_w),
      .surface_down_longwave_power_w =
          weighted(&GrayRadiationDiagnostics::surface_down_longwave_power_w),
      .surface_up_longwave_power_w =
          weighted(&GrayRadiationDiagnostics::surface_up_longwave_power_w),
      .atmospheric_shortwave_heating_power_w =
          weighted(&GrayRadiationDiagnostics::atmospheric_shortwave_heating_power_w),
      .atmospheric_longwave_heating_power_w =
          weighted(&GrayRadiationDiagnostics::atmospheric_longwave_heating_power_w),
      .surface_storage_rate_w =
          weighted(&GrayRadiationDiagnostics::surface_storage_rate_w),
      .sensible_to_atmosphere_power_w =
          weighted(&GrayRadiationDiagnostics::sensible_to_atmosphere_power_w),
      .internal_heat_power_w =
          weighted(&GrayRadiationDiagnostics::internal_heat_power_w),
      .interface_conservation_residual_w =
          weighted(&GrayRadiationDiagnostics::interface_conservation_residual_w),
      .dry_thermal_energy_rate_w =
          weighted(&GrayRadiationDiagnostics::dry_thermal_energy_rate_w),
      .rayleigh_drag_work_w = weighted(&GrayRadiationDiagnostics::rayleigh_drag_work_w),
  };
}

}  // namespace

DryHydrostaticDriver::DryHydrostaticDriver(ExperimentConfig c)
    : config_(std::move(c)),
      grid_(config_.grid.cells_per_panel, config_.planet.radius_m),
      coordinate_({config_.vertical.a_half_pa, config_.vertical.b_half},
                  config_.vertical.minimum_surface_pressure_pa,
                  config_.vertical.maximum_surface_pressure_pa,
                  config_.vertical.minimum_pressure_thickness_pa),
      orography_(make_surface_orography(config_.orography, grid_, config_.planet,
                                        config_.source_directory)) {
  if (config_.kind != ExperimentKind::kDryHydrostatic)
    throw std::invalid_argument("dry driver requires dry_hydrostatic config");
  if (config_.surface.has_value()) {
    surface_boundary_ =
        make_surface_boundary(*config_.surface, config_.orography, grid_,
                              config_.planet, config_.source_directory);
    orography_ = SurfaceOrography(
        std::vector<Real>(surface_boundary_->surface_geopotential_m2_s2().begin(),
                          surface_boundary_->surface_geopotential_m2_s2().end()),
        surface_boundary_->source_fingerprint());
  }
  // Reference-state well balancing is intentionally limited to the named hydrostatic
  // rest benchmark. Applying an initial-state correction to jets or forced cases would
  // silently modify their physical pressure force.
  if (config_.dry_hydrostatic.test_case == DryHydrostaticTestCase::kDcmip200Rest) {
    const auto reference_state = initial_state();
    pressure_reference_ = make_dry_hydrostatic_pressure_reference(
        grid_, diagnose(reference_state), config_.planet);
  }
  if (config_.dry_hydrostatic.time_integrator ==
      DryHydrostaticTimeIntegrator::kSemiImplicit) {
    if (!config_.semi_implicit.has_value())
      throw std::invalid_argument(
          "semi-implicit dry driver requires semi-implicit parameters");
    semi_implicit_reference_column_ = make_dry_hydrostatic_reference_column(
        coordinate_, config_.planet, *config_.semi_implicit);
    semi_implicit_fast_operator_ = make_dry_hydrostatic_fast_operator(
        coordinate_, config_.planet, *semi_implicit_reference_column_);
    semi_implicit_vertical_modes_ = make_dry_hydrostatic_vertical_modes(
        *semi_implicit_reference_column_, config_.planet,
        *semi_implicit_fast_operator_);
    semi_implicit_external_operator_ = make_dry_hydrostatic_external_mode_operator(
        *semi_implicit_reference_column_, *semi_implicit_vertical_modes_);
  }
}

void DryHydrostaticDriver::update_semi_implicit_reference(
    const DryHydrostaticState& state) const {
  // ADR 0016: keep the reference column horizontally uniform so a single vertical
  // structure eigendecomposition still serves every cell, but track the current state
  // so the reference-linear operator stays close to the true Jacobian.
  const auto levels = coordinate_.levels();
  const auto cells = grid_.cell_count();
  workspace_.column_potential_temperature.resize(levels);
  diagnose_dry_hydrostatic_state(
      state, coordinate_, config_.planet, orography_.surface_geopotential_m2_s2(),
      workspace_.derived, workspace_.vertical_geometry, workspace_.hydrostatic_column,
      workspace_.column_potential_temperature);
  std::vector<Real> profile(levels, 0.0);
  Real area = 0.0;
  Real surface_pressure = 0.0;
  for (std::size_t cell = 0; cell < cells; ++cell) {
    const Real cell_area = grid_.cells()[cell].area_m2;
    area += cell_area;
    surface_pressure += cell_area * state.surface_pressure_pa[cell];
    for (std::size_t level = 0; level < levels; ++level)
      profile[level] +=
          cell_area *
          workspace_.derived.temperature_k[dry_hydrostatic_offset(cell, level, levels)];
  }
  if (!(area > 0.0)) throw std::runtime_error("grid area is invalid");
  surface_pressure /= area;
  for (auto& value : profile) value /= area;
  semi_implicit_reference_column_ = make_dry_hydrostatic_reference_column(
      coordinate_, config_.planet, surface_pressure, profile);
  semi_implicit_fast_operator_ = make_dry_hydrostatic_fast_operator(
      coordinate_, config_.planet, *semi_implicit_reference_column_);
  semi_implicit_vertical_modes_ = make_dry_hydrostatic_vertical_modes(
      *semi_implicit_reference_column_, config_.planet, *semi_implicit_fast_operator_);
  semi_implicit_external_operator_ = make_dry_hydrostatic_external_mode_operator(
      *semi_implicit_reference_column_, *semi_implicit_vertical_modes_);
}
DryHydrostaticState DryHydrostaticDriver::initial_state() const {
  auto state =
      initialize_dry_hydrostatic_benchmark(config_, grid_, coordinate_, orography_);
  if (config_.physics.kind == PhysicsKind::kSurfaceEnergyBalance ||
      config_.physics.kind == PhysicsKind::kGrayRadiation)
    state.surface_temperature_k.assign(grid_.cell_count(),
                                       config_.surface->initial_temperature_k);
  return state;
}
DryHydrostaticDerived DryHydrostaticDriver::diagnose(
    const DryHydrostaticState& s) const {
  return diagnose_dry_hydrostatic_state(s, coordinate_, config_.planet,
                                        orography_.surface_geopotential_m2_s2());
}
DryHydrostaticSources DryHydrostaticDriver::diagnose_sources(
    const DryHydrostaticDerived& d) const {
  DryHydrostaticSources result;
  DryHydrostaticSourcesWorkspace workspace;
  if (pressure_reference_.has_value())
    dry_hydrostatic_sources(grid_, d, config_.planet, *pressure_reference_, result,
                            workspace);
  else
    dry_hydrostatic_sources(grid_, d, config_.planet, result, workspace);
  return result;
}
DryHydrostaticRhs DryHydrostaticDriver::rhs(const DryHydrostaticState& s) const {
  DryHydrostaticRhs result;
  rhs(s, result);
  return result;
}

void DryHydrostaticDriver::rhs(const DryHydrostaticState& s,
                               DryHydrostaticRhs& result) const {
  rhs_with_components(s, result, nullptr);
}

DryHydrostaticRhsComponents DryHydrostaticDriver::rhs_components(
    const DryHydrostaticState& s) const {
  DryHydrostaticRhs result;
  DryHydrostaticRhsComponents components;
  rhs_with_components(s, result, &components);
  return components;
}

void DryHydrostaticDriver::rhs_with_components(
    const DryHydrostaticState& s, DryHydrostaticRhs& result,
    DryHydrostaticRhsComponents* const components,
    const bool compute_fast_wave_cfl) const {
  auto& workspace = workspace_;
  workspace.column_potential_temperature.resize(coordinate_.levels());
  diagnose_dry_hydrostatic_state(
      s, coordinate_, config_.planet, orography_.surface_geopotential_m2_s2(),
      workspace.derived, workspace.vertical_geometry, workspace.hydrostatic_column,
      workspace.column_potential_temperature);
  const auto& d = workspace.derived;
  auto C = d.cells, K = d.levels;
  if (components != nullptr) resize_zero_rhs_components(*components, C, K);
  auto& h = workspace.horizontal_tendency;
  h.air_mass.assign(C * K, 0.0);
  h.momentum.assign(C * K, {});
  h.potential_temperature_mass.assign(C * K, 0.0);
  h.tracer_mass.assign(C * K, 0.0);
  // Per-cell Courant condition (ADR 0011): dt * sum_f(lambda_f * L_f) / A_cell <= cfl,
  // the same definition the transport and shallow-water solvers already use. The
  // previous per-edge form was about four times weaker on a quadrilateral cell.
  if (compute_fast_wave_cfl) workspace.face_fast_wave_speed_length.assign(C * K, 0.0);
  workspace.face_advective_speed_length.assign(C * K, 0.0);
  auto& face_fast_wave_speed_length = workspace.face_fast_wave_speed_length;
  auto& face_advective_speed_length = workspace.face_advective_speed_length;
  reconstruct_dry_hydrostatic_face_states(
      grid_, d, config_.dry_hydrostatic.reconstruction, config_.dry_hydrostatic.limiter,
      workspace.reconstruction, workspace.reconstruction_workspace,
      compute_fast_wave_cfl);
  const auto& reconstructed = workspace.reconstruction;
  for (const auto& e : grid_.edges()) {
    const auto& cached_edge = grid_.edge_cache()[e.id];
    auto l = cached_edge.left_cell;
    auto r = cached_edge.right_cell;
    const EdgeTangentBasis basis{cached_edge.normal, cached_edge.tangent};
    for (std::size_t k = 0; k < K; ++k) {
      const auto& face = reconstructed.at(e.id, k);
      auto f = rusanov_dry_hydrostatic_flux(
          face.left, face.right, basis, config_.planet.gas_constant_j_kg_k,
          config_.planet.heat_capacity_cp_j_kg_k, compute_fast_wave_cfl);
      auto add = [&](std::size_t c, Real sign) {
        auto n = dry_hydrostatic_offset(c, k, K);
        auto scale = sign * e.length_m / grid_.cells()[c].area_m2;
        h.air_mass[n] += scale * f.air_mass_kg_m_s;
        h.momentum[n] = h.momentum[n] + scale * f.momentum_kg_s2;
        h.potential_temperature_mass[n] +=
            scale * f.potential_temperature_mass_k_kg_m_s;
        h.tracer_mass[n] += scale * f.tracer_mass_kg_m_s;
        if (compute_fast_wave_cfl)
          face_fast_wave_speed_length[n] += e.length_m * f.maximum_wave_speed_m_s;
        face_advective_speed_length[n] += e.length_m * f.maximum_dissipation_speed_m_s;
      };
      add(l, -1);
      add(r, 1);
    }
  }
  Real fast_wave_dt = std::numeric_limits<Real>::infinity();
  Real advective_dt = std::numeric_limits<Real>::infinity();
  for (std::size_t c = 0; c < C; ++c) {
    for (std::size_t k = 0; k < K; ++k) {
      const auto n = dry_hydrostatic_offset(c, k, K);
      if (compute_fast_wave_cfl) {
        const Real fast_wave_denominator = face_fast_wave_speed_length[n];
        if (fast_wave_denominator > 0.0)
          fast_wave_dt = std::min(fast_wave_dt, config_.dry_hydrostatic.cfl *
                                                    grid_.cells()[c].area_m2 /
                                                    fast_wave_denominator);
      }
      const Real advective_denominator = face_advective_speed_length[n];
      if (advective_denominator > 0.0)
        advective_dt = std::min(advective_dt, config_.dry_hydrostatic.cfl *
                                                  grid_.cells()[c].area_m2 /
                                                  advective_denominator);
    }
  }
  couple_dry_hydrostatic_columns(
      s, d, h, coordinate_.coefficients().b_half, config_.planet.gravity_m_s2,
      config_.vertical.transport_scheme, config_.vertical.limiter,
      // Nonlinear minmod switches in theta transport degrade the cancellation between
      // thermodynamic transport and hydrostatic pressure work. Keep the configured
      // limiter for momentum and passive tracer, but use the smoother linear theta
      // reconstruction. Stage validation still rejects non-positive temperatures.
      VerticalLimiterKind::kNone, workspace.coupling, workspace.coupling_workspace);
  auto& coupled = workspace.coupling;
  if (components != nullptr) {
    components->horizontal_transport.tendency = h;
    components->vertical_transport.surface_pressure_pa_s =
        coupled.surface_pressure_pa_s;
    for (std::size_t n = 0; n < C * K; ++n) {
      components->vertical_transport.tendency.air_mass[n] =
          coupled.tendency.air_mass[n] - h.air_mass[n];
      components->vertical_transport.tendency.momentum[n] =
          coupled.tendency.momentum[n] - h.momentum[n];
      components->vertical_transport.tendency.potential_temperature_mass[n] =
          coupled.tendency.potential_temperature_mass[n] -
          h.potential_temperature_mass[n];
      components->vertical_transport.tendency.tracer_mass[n] =
          coupled.tendency.tracer_mass[n] - h.tracer_mass[n];
    }
  }
  Real vertical_dt = std::numeric_limits<Real>::infinity();
  for (std::size_t c = 0; c < C; ++c) {
    const auto begin = c * K;
    const auto interface_begin = c * (K + 1);
    const std::span<const Real> air_mass(d.air_mass_kg_m2.data() + begin, K);
    const std::span<const Real> interface_flux(
        coupled.interface_mass_flux_kg_m2_s.data() + interface_begin, K + 1);
    const std::span<const Real> horizontal_air_mass(h.air_mass.data() + begin, K);
    const Real unit_step_cfl =
        vertical_maximum_cfl(air_mass, interface_flux, horizontal_air_mass, 1.0);
    if (unit_step_cfl > 0.0)
      vertical_dt = std::min(vertical_dt, config_.vertical.cfl / unit_step_cfl);
  }
  if (pressure_reference_.has_value())
    dry_hydrostatic_sources(grid_, d, config_.planet, *pressure_reference_,
                            workspace.sources, workspace.sources_workspace);
  else
    dry_hydrostatic_sources(grid_, d, config_.planet, workspace.sources,
                            workspace.sources_workspace);
  const auto& sources = workspace.sources;
  if (components != nullptr) {
    components->pressure_gradient.tendency.momentum = sources.pressure_gradient_kg_m_s2;
    components->coriolis.tendency.momentum = sources.coriolis_kg_m_s2;
  }
  for (std::size_t n = 0; n < C * K; ++n)
    coupled.tendency.momentum[n] = coupled.tendency.momentum[n] +
                                   sources.pressure_gradient_kg_m_s2[n] +
                                   sources.coriolis_kg_m_s2[n];
  Real diffusion_rate = 0.0;
  Real diffusion_dt = std::numeric_limits<Real>::infinity();
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
    if (components != nullptr) {
      components->diffusion.tendency.momentum = diffusion.momentum;
      components->diffusion.tendency.potential_temperature_mass =
          diffusion.potential_temperature_mass;
      components->diffusion.tendency.tracer_mass = diffusion.tracer_mass;
    }
    diffusion_dt =
        stable_diffusion_time_step(grid_, config_.dry_hydrostatic.diffusion_kind,
                                   config_.dry_hydrostatic.diffusion_coefficient,
                                   std::numeric_limits<Real>::max());
  }
  HeldSuarezDiagnostics physics_diagnostics{};
  SurfaceEnergyDiagnostics surface_diagnostics{};
  GrayRadiationDiagnostics radiation_diagnostics{};
  Real surface_dt = std::numeric_limits<Real>::infinity();
  Real radiation_dt = std::numeric_limits<Real>::infinity();
  result.surface_temperature_k_s.clear();
  if (config_.physics.kind == PhysicsKind::kHeldSuarez) {
    held_suarez_tendency(grid_, coordinate_, d, s.surface_pressure_pa, config_.planet,
                         workspace.atmospheric_physics,
                         workspace.atmospheric_physics_workspace);
    const auto& physics = workspace.atmospheric_physics;
    if (components != nullptr) {
      components->physics.tendency.momentum = physics.horizontal_momentum_mass_kg_m_s2;
      components->physics.tendency.potential_temperature_mass =
          physics.potential_temperature_mass_k_kg_m2_s;
    }
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
    surface_energy_tendency(grid_, *surface_boundary_, s.surface_temperature_k, d,
                            s.surface_pressure_pa, config_.planet, *config_.surface,
                            orbit_state, workspace.surface_physics);
    const auto& physics = workspace.surface_physics;
    if (components != nullptr) {
      components->physics.tendency.momentum = physics.horizontal_momentum_mass_kg_m_s2;
      components->physics.tendency.potential_temperature_mass =
          physics.potential_temperature_mass_k_kg_m2_s;
      components->physics.surface_temperature_k_s = physics.surface_temperature_k_s;
    }
    for (std::size_t n = 0; n < C * K; ++n) {
      coupled.tendency.momentum[n] =
          coupled.tendency.momentum[n] + physics.horizontal_momentum_mass_kg_m_s2[n];
      coupled.tendency.potential_temperature_mass[n] +=
          physics.potential_temperature_mass_k_kg_m2_s[n];
    }
    result.surface_temperature_k_s = physics.surface_temperature_k_s;
    surface_diagnostics = physics.diagnostics;
    surface_dt = physics.stable_time_step_s;
  } else if (config_.physics.kind == PhysicsKind::kPlanetaryNewtonian) {
    std::optional<OrbitState> orbit_state;
    if (config_.physics.geometry == ForcingGeometry::kSubstellar)
      orbit_state = evaluate_orbit(*config_.orbit, config_.planet.rotation_rate_rad_s,
                                   s.time_s - config_.run.start_time_s);
    planetary_newtonian_tendency(
        grid_, coordinate_, d, s.surface_pressure_pa, config_.planet,
        config_.physics.geometry, orbit_state ? &*orbit_state : nullptr,
        workspace.atmospheric_physics, workspace.atmospheric_physics_workspace);
    const auto& physics = workspace.atmospheric_physics;
    if (components != nullptr) {
      components->physics.tendency.momentum = physics.horizontal_momentum_mass_kg_m_s2;
      components->physics.tendency.potential_temperature_mass =
          physics.potential_temperature_mass_k_kg_m2_s;
    }
    for (std::size_t n = 0; n < C * K; ++n) {
      coupled.tendency.momentum[n] =
          coupled.tendency.momentum[n] + physics.horizontal_momentum_mass_kg_m_s2[n];
      coupled.tendency.potential_temperature_mass[n] +=
          physics.potential_temperature_mass_k_kg_m2_s[n];
    }
    physics_diagnostics = physics.diagnostics;
  } else if (config_.physics.kind == PhysicsKind::kGrayRadiation) {
    const auto orbit_state =
        evaluate_orbit(*config_.orbit, config_.planet.rotation_rate_rad_s,
                       s.time_s - config_.run.start_time_s);
    gray_radiation_tendency(grid_, coordinate_, *surface_boundary_,
                            s.surface_temperature_k, d, s.surface_pressure_pa,
                            config_.planet, *config_.surface, *config_.radiation,
                            orbit_state, workspace.gray_radiation_physics,
                            workspace.gray_radiation_workspace);
    const auto& physics = workspace.gray_radiation_physics;
    if (components != nullptr) {
      components->physics.tendency.momentum = physics.horizontal_momentum_mass_kg_m_s2;
      components->physics.tendency.potential_temperature_mass =
          physics.potential_temperature_mass_k_kg_m2_s;
      components->physics.surface_temperature_k_s = physics.surface_temperature_k_s;
    }
    for (std::size_t n = 0; n < C * K; ++n) {
      coupled.tendency.momentum[n] =
          coupled.tendency.momentum[n] + physics.horizontal_momentum_mass_kg_m_s2[n];
      coupled.tendency.potential_temperature_mass[n] +=
          physics.potential_temperature_mass_k_kg_m2_s[n];
    }
    result.surface_temperature_k_s = physics.surface_temperature_k_s;
    radiation_diagnostics = physics.diagnostics;
    radiation_dt = physics.stable_time_step_s;
    physics_diagnostics.thermal_energy_rate_w =
        physics.diagnostics.dry_thermal_energy_rate_w;
    physics_diagnostics.rayleigh_drag_work_w = physics.diagnostics.rayleigh_drag_work_w;
  }
  result.surface_pressure_pa_s = coupled.surface_pressure_pa_s;
  result.tendency = coupled.tendency;
  result.horizontal_fast_wave_stable_time_step_s = fast_wave_dt;
  result.horizontal_advective_stable_time_step_s = advective_dt;
  result.diffusion_stable_time_step_s = diffusion_dt;
  result.horizontal_stable_time_step_s = std::min(fast_wave_dt, diffusion_dt);
  result.vertical_stable_time_step_s = vertical_dt;
  result.surface_stable_time_step_s = surface_dt;
  result.radiation_stable_time_step_s = radiation_dt;
  result.maximum_continuity_residual_pa_s = coupled.maximum_continuity_residual_pa_s;
  result.physics_diagnostics = physics_diagnostics;
  result.surface_diagnostics = surface_diagnostics;
  result.radiation_diagnostics = radiation_diagnostics;
  result.diffusion_kinetic_energy_rate_w = diffusion_rate;
}
void DryHydrostaticDriver::advance(DryHydrostaticState& s, const Real end,
                                   const DryHydrostaticObserver& obs,
                                   const DryHydrostaticCancel& cancel) const {
  const std::uint64_t initial_observer_step = s.step;
  const auto observe = [&](const DryHydrostaticState& state,
                           const DryHydrostaticStepDiagnostics& step) {
    if (!obs) return;
    const bool needs_sample = state.step == initial_observer_step ||
                              state.step % config_.diagnostics.interval_steps == 0 ||
                              state.time_s >= end;
    if (!needs_sample) {
      obs(state, nullptr, step);
      return;
    }
    workspace_.column_potential_temperature.resize(coordinate_.levels());
    diagnose_dry_hydrostatic_state(
        state, coordinate_, config_.planet, orography_.surface_geopotential_m2_s2(),
        workspace_.derived, workspace_.vertical_geometry, workspace_.hydrostatic_column,
        workspace_.column_potential_temperature);
    const auto& derived = workspace_.derived;
    auto sampled = step;
    if (config_.physics.kind == PhysicsKind::kHeldSuarez) {
      held_suarez_tendency(grid_, coordinate_, derived, state.surface_pressure_pa,
                           config_.planet, workspace_.atmospheric_physics,
                           workspace_.atmospheric_physics_workspace);
      sampled.physics_rates = workspace_.atmospheric_physics.diagnostics;
    } else if (config_.physics.kind == PhysicsKind::kPlanetaryNewtonian) {
      std::optional<OrbitState> orbit_state;
      if (config_.physics.geometry == ForcingGeometry::kSubstellar)
        orbit_state = evaluate_orbit(*config_.orbit, config_.planet.rotation_rate_rad_s,
                                     state.time_s - config_.run.start_time_s);
      planetary_newtonian_tendency(
          grid_, coordinate_, derived, state.surface_pressure_pa, config_.planet,
          config_.physics.geometry, orbit_state ? &*orbit_state : nullptr,
          workspace_.atmospheric_physics, workspace_.atmospheric_physics_workspace);
      sampled.physics_rates = workspace_.atmospheric_physics.diagnostics;
    } else if (config_.physics.kind == PhysicsKind::kSurfaceEnergyBalance) {
      const auto orbit_state =
          evaluate_orbit(*config_.orbit, config_.planet.rotation_rate_rad_s,
                         state.time_s - config_.run.start_time_s);
      surface_energy_tendency(grid_, *surface_boundary_, state.surface_temperature_k,
                              derived, state.surface_pressure_pa, config_.planet,
                              *config_.surface, orbit_state,
                              workspace_.surface_physics);
      sampled.surface_rates = workspace_.surface_physics.diagnostics;
    } else if (config_.physics.kind == PhysicsKind::kGrayRadiation) {
      const auto orbit_state =
          evaluate_orbit(*config_.orbit, config_.planet.rotation_rate_rad_s,
                         state.time_s - config_.run.start_time_s);
      gray_radiation_tendency(
          grid_, coordinate_, *surface_boundary_, state.surface_temperature_k, derived,
          state.surface_pressure_pa, config_.planet, *config_.surface,
          *config_.radiation, orbit_state, workspace_.gray_radiation_physics,
          workspace_.gray_radiation_workspace);
      sampled.radiation_rates = workspace_.gray_radiation_physics.diagnostics;
    }
    obs(state, &derived, sampled);
  };
  observe(s, {});
  std::vector<Vec3> centres;
  centres.reserve(grid_.cell_count());
  for (const auto& cell : grid_.cells()) centres.push_back(cell.center);
  const auto project_momentum = [&](DryHydrostaticState& stage) {
    for (std::size_t n = 0; n < stage.horizontal_momentum_mass_kg_m_s.size(); ++n) {
      const auto cell = n / coordinate_.levels();
      stage.horizontal_momentum_mass_kg_m_s[n] = project_tangent(
          stage.horizontal_momentum_mass_kg_m_s[n], grid_.cells()[cell].center);
    }
  };
  const auto diagnose_and_validate = [&](const DryHydrostaticState& stage) {
    workspace_.column_potential_temperature.resize(coordinate_.levels());
    diagnose_dry_hydrostatic_state(
        stage, coordinate_, config_.planet, orography_.surface_geopotential_m2_s2(),
        workspace_.derived, workspace_.vertical_geometry, workspace_.hydrostatic_column,
        workspace_.column_potential_temperature);
    validate_dry_hydrostatic_state(stage, workspace_.derived, centres,
                                   config_.vertical.minimum_surface_pressure_pa,
                                   config_.vertical.maximum_surface_pressure_pa,
                                   config_.vertical.temperature_floor_k);
  };
  if (config_.dry_hydrostatic.time_integrator ==
      DryHydrostaticTimeIntegrator::kSemiImplicit) {
    const auto& parameters = *config_.semi_implicit;
    const bool reference_is_per_step =
        parameters.reference_update == SemiImplicitReferenceUpdate::kPerStep;
    const GmresOptions linear_options{
        .restart = static_cast<std::size_t>(parameters.gmres_restart),
        .maximum_iterations =
            static_cast<std::size_t>(parameters.linear_maximum_iterations),
        .relative_tolerance = parameters.linear_relative_tolerance,
        .absolute_tolerance = parameters.linear_absolute_tolerance};
    DryHydrostaticRhs initial_rhs;
    DryHydrostaticRhs candidate_rhs;
    bool initial_rhs_is_current = false;
    const bool rhs_is_time_independent =
        config_.physics.kind != PhysicsKind::kPlanetaryNewtonian &&
        config_.physics.kind != PhysicsKind::kSurfaceEnergyBalance &&
        config_.physics.kind != PhysicsKind::kGrayRadiation;
    while (s.time_s < end) {
      if (cancel && cancel()) return;
      const auto step_wall_start = std::chrono::steady_clock::now();
      Real rhs_wall_seconds = 0.0;
      Real linear_solve_wall_seconds = 0.0;
      std::size_t linear_iterations_total = 0;
      std::size_t linear_iterations_maximum = 0;
      Real linear_relative_residual_maximum = 0.0;
      std::size_t accepted_selected_modes = 0;
      std::size_t accepted_nonlinear_iterations = 0;
      Real accepted_nonlinear_residual = 0.0;
      std::size_t retry_count = 0;
      const DryHydrostaticState initial = s;
      if (reference_is_per_step) update_semi_implicit_reference(initial);
      const auto& reference = *semi_implicit_reference_column_;
      const auto& fast_operator = *semi_implicit_fast_operator_;
      const auto& vertical_modes = *semi_implicit_vertical_modes_;
      const auto& external_operator = *semi_implicit_external_operator_;
      const Real requested_dt = std::min(config_.run.time_step_s, end - initial.time_s);
      const auto timed_rhs = [&](const DryHydrostaticState& state,
                                 DryHydrostaticRhs& result) {
        const auto start = std::chrono::steady_clock::now();
        rhs_with_components(state, result, nullptr, false);
        rhs_wall_seconds +=
            std::chrono::duration<Real>(std::chrono::steady_clock::now() - start)
                .count();
      };
      if (!initial_rhs_is_current) timed_rhs(initial, initial_rhs);
      Real dt = std::min({config_.run.time_step_s,
                          semi_implicit_explicit_limit(config_, initial_rhs),
                          end - initial.time_s});
      dt = pressure_limited_time_step(config_, initial, initial_rhs, dt);
      if (!(dt > 0.0) || !std::isfinite(dt))
        throw std::runtime_error("semi-implicit dry time step is invalid");
      if (dt < parameters.minimum_time_step_s && dt < end - initial.time_s)
        throw std::runtime_error(
            "semi-implicit dry time step is below the configured minimum");
      const Real minimum_retry_time_step_s =
          std::min(parameters.minimum_time_step_s, end - initial.time_s);
      const auto attempt_step = [&](const Real attempted_dt) {
        const auto selected_modes = select_implicit_vertical_modes(
            grid_, vertical_modes, attempted_dt, parameters.wave_cfl_threshold,
            static_cast<std::size_t>(parameters.maximum_implicit_modes));
        accepted_selected_modes = selected_modes.size();
        DryHydrostaticState candidate = initial;
        candidate.time_s = initial.time_s + attempted_dt;
        // ADR 0016 / Benard (2003): an ICI scheme performs a fixed number of
        // quasi-Newton iterations. The residual is a diagnostic, never an acceptance
        // test; only a non-finite or growing residual rejects the attempt.
        Real initial_nonlinear_residual = 0.0;
        accepted_nonlinear_iterations =
            static_cast<std::size_t>(parameters.nonlinear_iterations);
        for (Index iteration = 0; iteration < parameters.nonlinear_iterations;
             ++iteration) {
          const DryHydrostaticRhs* rhs_at_candidate = nullptr;
          if (iteration == 0 && rhs_is_time_independent) {
            rhs_at_candidate = &initial_rhs;
          } else {
            timed_rhs(candidate, candidate_rhs);
            rhs_at_candidate = &candidate_rhs;
          }
          const auto residual = crank_nicolson_residual(
              initial, candidate, initial_rhs, *rhs_at_candidate, attempted_dt,
              parameters.implicit_weight, centres, coordinate_.levels());
          const auto scaled_residual = scaled_crank_nicolson_residual(
              residual, initial, reference, external_operator);
          if (!std::isfinite(scaled_residual.norm))
            throw std::runtime_error("Crank-Nicolson residual is not finite");
          if (iteration == 0) initial_nonlinear_residual = scaled_residual.norm;
          accepted_nonlinear_residual = scaled_residual.norm;
          const auto correction_rhs = negative_fast_residual(residual);
          const auto solve_start = std::chrono::steady_clock::now();
          const auto solve = solve_dry_hydrostatic_modal_correction(
              grid_, config_.planet, fast_operator, vertical_modes, selected_modes,
              parameters.implicit_weight * attempted_dt, correction_rhs, linear_options,
              semi_implicit_workspace_.correction, semi_implicit_workspace_);
          linear_solve_wall_seconds +=
              std::chrono::duration<Real>(std::chrono::steady_clock::now() -
                                          solve_start)
                  .count();
          linear_iterations_total += solve.linear_iterations_total;
          linear_iterations_maximum =
              std::max(linear_iterations_maximum, solve.linear_iterations_maximum);
          linear_relative_residual_maximum = std::max(
              linear_relative_residual_maximum, solve.linear_relative_residual_maximum);
          if (!solve.all_converged || solve.equation_residual_norm >
                                          10.0 * parameters.linear_relative_tolerance)
            throw std::runtime_error(
                modal_failure_message(solve.equation_residual_norm));
          add_semi_implicit_correction(candidate, semi_implicit_workspace_.correction,
                                       residual);
          project_momentum(candidate);
          if (!pressure_is_in_range(config_, candidate) ||
              !surface_temperature_is_positive(candidate))
            throw std::runtime_error("correction violates a prognostic invariant");
          diagnose_and_validate(candidate);
        }
        // Final evaluation at the accepted candidate: it supplies the explicit CFL
        // check, the step energy attribution, the reported residual, and the
        // first-same- as-last tendency reused as the next step's `initial_rhs`.
        timed_rhs(candidate, candidate_rhs);
        if (attempted_dt > semi_implicit_explicit_limit(config_, candidate_rhs))
          throw std::runtime_error("candidate violates an explicit CFL constraint");
        const auto final_residual = crank_nicolson_residual(
            initial, candidate, initial_rhs, candidate_rhs, attempted_dt,
            parameters.implicit_weight, centres, coordinate_.levels());
        const auto final_scaled = scaled_crank_nicolson_residual(
            final_residual, initial, reference, external_operator);
        accepted_nonlinear_residual = final_scaled.norm;
        if (!std::isfinite(final_scaled.norm) ||
            final_scaled.norm > initial_nonlinear_residual)
          throw std::runtime_error(nonlinear_failure_message(
              initial_nonlinear_residual, final_scaled.norm, final_scaled.rms_norm,
              final_scaled.component, attempted_dt));
        candidate.step = initial.step + 1;
        diagnose_and_validate(candidate);
        return candidate;
      };

      DryHydrostaticState candidate;
      while (true) {
        try {
          candidate = attempt_step(dt);
          break;
        } catch (const std::runtime_error& error) {
          const Real retry_time_step_s = 0.5 * dt;
          if (!(retry_time_step_s > 0.0) ||
              initial.time_s + retry_time_step_s == initial.time_s ||
              retry_time_step_s < minimum_retry_time_step_s) {
            throw std::runtime_error(
                "semi-implicit dry step failed at the configured minimum: " +
                std::string(error.what()));
          }
          dt = retry_time_step_s;
          ++retry_count;
        }
      }

      s = std::move(candidate);
      const Real explicit_weight = 1.0 - parameters.implicit_weight;
      DryHydrostaticStepDiagnostics step{
          .thermal_energy_contribution_j =
              dt *
              (explicit_weight * initial_rhs.physics_diagnostics.thermal_energy_rate_w +
               parameters.implicit_weight *
                   candidate_rhs.physics_diagnostics.thermal_energy_rate_w),
          .rayleigh_drag_energy_contribution_j =
              dt *
              (explicit_weight * initial_rhs.physics_diagnostics.rayleigh_drag_work_w +
               parameters.implicit_weight *
                   candidate_rhs.physics_diagnostics.rayleigh_drag_work_w),
          .diffusion_energy_contribution_j =
              dt * (explicit_weight * initial_rhs.diffusion_kinetic_energy_rate_w +
                    parameters.implicit_weight *
                        candidate_rhs.diffusion_kinetic_energy_rate_w),
          .requested_time_step_s = requested_dt,
          .accepted_time_step_s = dt,
          .advective_cfl =
              std::max(courant_from_stable_time_step(
                           dt, config_.dry_hydrostatic.cfl,
                           initial_rhs.horizontal_advective_stable_time_step_s),
                       courant_from_stable_time_step(
                           dt, config_.dry_hydrostatic.cfl,
                           candidate_rhs.horizontal_advective_stable_time_step_s)),
          .implicit_wave_courant =
              maximum_vertical_mode_courant(grid_, vertical_modes, dt),
          .vertical_cfl = std::max(
              courant_from_stable_time_step(dt, config_.vertical.cfl,
                                            initial_rhs.vertical_stable_time_step_s),
              courant_from_stable_time_step(dt, config_.vertical.cfl,
                                            candidate_rhs.vertical_stable_time_step_s)),
          .selected_implicit_modes = accepted_selected_modes,
          .linear_iterations_total = linear_iterations_total,
          .linear_iterations_maximum = linear_iterations_maximum,
          .linear_relative_residual_maximum = linear_relative_residual_maximum,
          .nonlinear_iterations = accepted_nonlinear_iterations,
          .nonlinear_relative_residual = accepted_nonlinear_residual,
          .retry_count = retry_count,
          .wall_seconds_rhs = rhs_wall_seconds,
          .wall_seconds_linear_solve = linear_solve_wall_seconds};
      if (config_.physics.kind == PhysicsKind::kGrayRadiation)
        step.radiation_rates = weighted_radiation_diagnostics(
            initial_rhs.radiation_diagnostics, explicit_weight,
            candidate_rhs.radiation_diagnostics, parameters.implicit_weight);
      if (config_.physics.kind == PhysicsKind::kSurfaceEnergyBalance) {
        const auto weighted = [&](const Real first, const Real second) {
          return explicit_weight * first + parameters.implicit_weight * second;
        };
        const SurfaceEnergyDiagnostics rates{
            .absorbed_stellar_power_w =
                weighted(initial_rhs.surface_diagnostics.absorbed_stellar_power_w,
                         candidate_rhs.surface_diagnostics.absorbed_stellar_power_w),
            .internal_heat_power_w =
                weighted(initial_rhs.surface_diagnostics.internal_heat_power_w,
                         candidate_rhs.surface_diagnostics.internal_heat_power_w),
            .outgoing_longwave_power_w =
                weighted(initial_rhs.surface_diagnostics.outgoing_longwave_power_w,
                         candidate_rhs.surface_diagnostics.outgoing_longwave_power_w),
            .sensible_to_atmosphere_power_w = weighted(
                initial_rhs.surface_diagnostics.sensible_to_atmosphere_power_w,
                candidate_rhs.surface_diagnostics.sensible_to_atmosphere_power_w),
            .surface_storage_rate_w =
                weighted(initial_rhs.surface_diagnostics.surface_storage_rate_w,
                         candidate_rhs.surface_diagnostics.surface_storage_rate_w)};
        step.surface_budget = integrate_surface_energy_budget(
            grid_, *surface_boundary_, *config_.surface, initial.surface_temperature_k,
            s.surface_temperature_k, dt, rates, rates, rates);
      }
      std::swap(initial_rhs, candidate_rhs);
      initial_rhs_is_current = true;
      step.wall_seconds_total = std::chrono::duration<Real>(
                                    std::chrono::steady_clock::now() - step_wall_start)
                                    .count();
      observe(s, step);
    }
    return;
  }
  DryHydrostaticRhs rhs1;
  DryHydrostaticRhs rhs2;
  DryHydrostaticRhs rhs3;
  while (s.time_s < end) {
    if (cancel && cancel()) return;
    const DryHydrostaticState initial = s;
    rhs(initial, rhs1);
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
      project_momentum(stage1);
      rhs(stage1, rhs2);
      validate_dry_hydrostatic_state(stage1, workspace_.derived, centres,
                                     config_.vertical.minimum_surface_pressure_pa,
                                     config_.vertical.maximum_surface_pressure_pa,
                                     config_.vertical.temperature_floor_k);
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
      project_momentum(stage2);
      rhs(stage2, rhs3);
      validate_dry_hydrostatic_state(stage2, workspace_.derived, centres,
                                     config_.vertical.minimum_surface_pressure_pa,
                                     config_.vertical.maximum_surface_pressure_pa,
                                     config_.vertical.temperature_floor_k);
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
      project_momentum(next);
      diagnose_and_validate(next);
      s = std::move(next);
      DryHydrostaticStepDiagnostics step{
          .thermal_energy_contribution_j =
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
                    2.0 * rhs3.diffusion_kinetic_energy_rate_w / 3.0)};
      if (config_.physics.kind == PhysicsKind::kGrayRadiation)
        step.radiation_rates = weighted_radiation_diagnostics(
            rhs1.radiation_diagnostics, 1.0 / 6.0, rhs2.radiation_diagnostics,
            1.0 / 6.0, rhs3.radiation_diagnostics, 2.0 / 3.0);
      if (config_.physics.kind == PhysicsKind::kSurfaceEnergyBalance) {
        step.surface_budget = integrate_surface_energy_budget(
            grid_, *surface_boundary_, *config_.surface, initial.surface_temperature_k,
            s.surface_temperature_k, dt, rhs1.surface_diagnostics,
            rhs2.surface_diagnostics, rhs3.surface_diagnostics);
      }
      observe(s, step);
      break;
    }
  }
}
}  // namespace mps
