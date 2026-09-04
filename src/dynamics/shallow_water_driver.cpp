#include "myplanetsim/dynamics/shallow_water_driver.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "myplanetsim/core/validation.hpp"
#include "myplanetsim/dynamics/shallow_water_benchmarks.hpp"
#include "myplanetsim/dynamics/shallow_water_compatible.hpp"
#include "myplanetsim/dynamics/shallow_water_diffusion.hpp"
#include "myplanetsim/dynamics/shallow_water_initial_edits.hpp"
#include "myplanetsim/dynamics/shallow_water_rhs.hpp"
#include "myplanetsim/dynamics/surface_orography.hpp"

namespace mps {
namespace {

[[nodiscard]] Vec3 rotation_vector(const ExperimentConfig& config) {
  if (config.shallow_water.test_case == ShallowWaterTestCase::kWilliamson2 ||
      config.shallow_water.test_case == ShallowWaterTestCase::kWilliamson5 ||
      config.shallow_water.test_case == ShallowWaterTestCase::kWilliamson6 ||
      config.shallow_water.test_case == ShallowWaterTestCase::kGalewsky) {
    return config.planet.rotation_rate_rad_s *
           normalize(Vec3{config.shallow_water.flow_axis_x,
                          config.shallow_water.flow_axis_y,
                          config.shallow_water.flow_axis_z});
  }
  return {0.0, 0.0, config.planet.rotation_rate_rad_s};
}

void project_and_validate_stage(const CubedSphereGrid& grid, ShallowWaterState& state,
                                const Real depth_floor_m, const int stage,
                                const Real cfl) {
  project_shallow_water_momentum(grid, state);
  try {
    validate_shallow_water_state(grid, state, depth_floor_m);
  } catch (const std::exception& error) {
    const auto minimum = std::ranges::min_element(state.depth);
    std::ostringstream message;
    message << "shallow-water stage " << stage << " failed at CFL " << cfl
            << " with minimum depth " << *minimum << ": " << error.what();
    throw std::runtime_error(message.str());
  }
}

[[nodiscard]] ShallowWaterState euler_update(const ShallowWaterState& base,
                                             const ShallowWaterTendency& tendency,
                                             const Real scale) {
  ShallowWaterState result = base;
  for (std::size_t cell = 0; cell < result.cell_count(); ++cell) {
    result.depth[cell] += scale * tendency.depth[cell];
    result.momentum[cell] = result.momentum[cell] + scale * tendency.momentum[cell];
  }
  return result;
}

[[nodiscard]] ShallowWaterState convex_update(const ShallowWaterState& initial,
                                              const Real initial_weight,
                                              const ShallowWaterState& stage,
                                              const ShallowWaterTendency& tendency,
                                              const Real stage_weight,
                                              const Real time_step_s) {
  ShallowWaterState result = initial;
  for (std::size_t cell = 0; cell < result.cell_count(); ++cell) {
    result.depth[cell] =
        initial_weight * initial.depth[cell] +
        stage_weight * (stage.depth[cell] + time_step_s * tendency.depth[cell]);
    result.momentum[cell] =
        initial_weight * initial.momentum[cell] +
        stage_weight * (stage.momentum[cell] + time_step_s * tendency.momentum[cell]);
  }
  return result;
}

void ssp_rk3_step(const CubedSphereGrid& grid, ShallowWaterState& state,
                  const Real time_step_s, const ShallowWaterParameters& parameters,
                  const CubedSphereDualTopology* dual, const Real gravity_m_s2,
                  const Vec3 rotation_vector_rad_s, const Real cfl,
                  const std::span<const Real> surface_geopotential_m2_s2) {
  const auto assemble = [&](const ShallowWaterState& stage) {
    if (dual != nullptr) {
      return assemble_compatible_shallow_water_rhs(grid, *dual, stage, parameters,
                                                   gravity_m_s2, rotation_vector_rad_s);
    }
    return assemble_shallow_water_rhs(grid, stage, parameters, gravity_m_s2,
                                      rotation_vector_rad_s,
                                      surface_geopotential_m2_s2);
  };
  const ShallowWaterState initial = state;
  const auto rhs1 = assemble(initial);
  auto stage1 = euler_update(initial, rhs1.total, time_step_s);
  project_and_validate_stage(grid, stage1, parameters.depth_floor_m, 1, cfl);

  const auto rhs2 = assemble(stage1);
  auto stage2 = convex_update(initial, 0.75, stage1, rhs2.total, 0.25, time_step_s);
  project_and_validate_stage(grid, stage2, parameters.depth_floor_m, 2, cfl);

  const auto rhs3 = assemble(stage2);
  state = convex_update(initial, 1.0 / 3.0, stage2, rhs3.total, 2.0 / 3.0, time_step_s);
  project_and_validate_stage(grid, state, parameters.depth_floor_m, 3, cfl);
}

}  // namespace

Real stable_shallow_water_time_step(const CubedSphereGrid& grid,
                                    const ShallowWaterState& state,
                                    const Real gravity_m_s2, const Real cfl,
                                    const Real maximum_time_step_s) {
  require_positive(gravity_m_s2, "shallow-water gravity");
  if (!(cfl > 0.0 && cfl <= 1.0) || !std::isfinite(cfl)) {
    throw std::invalid_argument("shallow-water CFL must be in (0, 1]");
  }
  require_positive(maximum_time_step_s, "maximum shallow-water time step");
  validate_shallow_water_state(grid, state, 0.0);
  Real time_step = maximum_time_step_s;
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
    Real denominator = 0.0;
    for (const std::size_t edge_id : grid.cell_edges(grid.cell_id(cell))) {
      const auto& edge = grid.edge(edge_id);
      const Vec3 normal = normalize(project_tangent(
          static_cast<Real>(grid.edge_sign_for_cell(edge_id, grid.cell_id(cell))) *
              edge.outward_normal_from_left,
          grid.cells()[cell].center));
      denominator += edge.length_m * (std::abs(dot(state.velocity(cell), normal)) +
                                      std::sqrt(gravity_m_s2 * state.depth[cell]));
    }
    time_step = std::min(time_step, cfl * grid.cells()[cell].area_m2 / denominator);
  }
  return time_step;
}

Real shallow_water_cfl_number(const CubedSphereGrid& grid,
                              const ShallowWaterState& state, const Real gravity_m_s2,
                              const Real time_step_s) {
  require_positive(gravity_m_s2, "shallow-water gravity");
  require_positive(time_step_s, "shallow-water time step");
  validate_shallow_water_state(grid, state, 0.0);
  Real maximum_cfl = 0.0;
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
    Real denominator = 0.0;
    for (const std::size_t edge_id : grid.cell_edges(grid.cell_id(cell))) {
      const auto& edge = grid.edge(edge_id);
      const Vec3 normal = normalize(project_tangent(
          static_cast<Real>(grid.edge_sign_for_cell(edge_id, grid.cell_id(cell))) *
              edge.outward_normal_from_left,
          grid.cells()[cell].center));
      denominator += edge.length_m * (std::abs(dot(state.velocity(cell), normal)) +
                                      std::sqrt(gravity_m_s2 * state.depth[cell]));
    }
    maximum_cfl =
        std::max(maximum_cfl, time_step_s * denominator / grid.cells()[cell].area_m2);
  }
  return maximum_cfl;
}

ShallowWaterResult run_shallow_water(
    const ExperimentConfig& config, std::optional<ShallowWaterState> initial_state,
    const std::optional<std::uint64_t> stop_after_step,
    const std::span<const InitialConditionEditV1> initial_edits,
    const ShallowWaterRunHooks hooks) {
  config.validate();
  if (config.kind != ExperimentKind::kShallowWater) {
    throw std::invalid_argument("shallow-water driver requires shallow_water config");
  }
  const CubedSphereGrid grid(config.grid.cells_per_panel, config.planet.radius_m);
  const auto orography = make_surface_orography(
      config.orography, grid, config.planet.gravity_m_s2);
  const Vec3 omega = rotation_vector(config);
  std::optional<CubedSphereDualTopology> dual;
  if (config.shallow_water.scheme == ShallowWaterScheme::kCompatible) {
    dual.emplace(grid);
  }
  const ShallowWaterState reference = make_shallow_water_initial_state(grid, config);
  ShallowWaterState state =
      initial_state.has_value() ? std::move(*initial_state) : reference;
  if (!initial_edits.empty()) {
    if (initial_state.has_value()) {
      throw std::invalid_argument("initial edits cannot be applied to a restart state");
    }
    const auto edit_diagnostics = apply_initial_condition_edits(
        grid, state, initial_edits, config.shallow_water.depth_floor_m);
    (void)edit_diagnostics;
  }
  validate_shallow_water_state(grid, state, config.shallow_water.depth_floor_m);
  const auto initial_diagnostics = diagnostics::diagnose_shallow_water(
      grid, initial_edits.empty() ? reference : state, config.planet.gravity_m_s2,
      omega);
  if (state.time_s < config.run.start_time_s || state.time_s > config.run.end_time_s) {
    throw std::invalid_argument("shallow-water restart time is outside configured run");
  }
  Real maximum_cfl = 0.0;
  std::vector<ShallowWaterSample> samples;
  const auto sample = [&](const ShallowWaterState& sampled) {
    samples.push_back({.time_s = sampled.time_s,
                       .step = sampled.step,
                       .invariants = diagnostics::diagnose_shallow_water(
                           grid, sampled, config.planet.gravity_m_s2, omega)});
  };
  std::uint64_t last_notified_frame_step = std::numeric_limits<std::uint64_t>::max();
  const auto notify_frame = [&](const ShallowWaterState& sampled) {
    if (hooks.on_frame != nullptr) {
      hooks.on_frame(sampled, hooks.observer_context);
      last_notified_frame_step = sampled.step;
    }
  };
  sample(state);
  notify_frame(state);
  while (state.time_s < config.run.end_time_s &&
         (!stop_after_step.has_value() || state.step < *stop_after_step)) {
    if (hooks.is_cancelled != nullptr &&
        hooks.is_cancelled(hooks.cancellation_context)) {
      break;
    }
    Real time_step = stable_shallow_water_time_step(
        grid, state, config.planet.gravity_m_s2, config.shallow_water.cfl,
        config.run.time_step_s);
    time_step = stable_diffusion_time_step(grid, config.shallow_water.diffusion_kind,
                                           config.shallow_water.diffusion_coefficient,
                                           time_step);
    time_step = std::min(time_step, config.run.end_time_s - state.time_s);
    const Real actual_cfl =
        shallow_water_cfl_number(grid, state, config.planet.gravity_m_s2, time_step);
    maximum_cfl = std::max(maximum_cfl, actual_cfl);
    ssp_rk3_step(grid, state, time_step, config.shallow_water,
                 dual.has_value() ? &*dual : nullptr, config.planet.gravity_m_s2, omega,
                 actual_cfl, orography.surface_geopotential_m2_s2());
    state.time_s += time_step;
    ++state.step;
    if (state.step % config.diagnostics.interval_steps == 0) {
      sample(state);
      notify_frame(state);
    }
  }
  if (samples.back().step != state.step) {
    sample(state);
  }
  if (hooks.on_frame != nullptr) {
    if (last_notified_frame_step != state.step) {
      notify_frame(state);
    }
  }
  const auto final_diagnostics = diagnostics::diagnose_shallow_water(
      grid, state, config.planet.gravity_m_s2, omega);
  const bool reached_end_time = state.time_s == config.run.end_time_s;
  return {.state = std::move(state),
          .initial_diagnostics = initial_diagnostics,
          .final_diagnostics = final_diagnostics,
          .reached_end_time = reached_end_time,
          .maximum_cfl = maximum_cfl,
          .samples = std::move(samples)};
}

}  // namespace mps
