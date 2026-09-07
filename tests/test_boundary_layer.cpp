#include <algorithm>
#include <array>
#include <cmath>
#include <numbers>
#include <utility>
#include <vector>

#include "myplanetsim/physics/boundary_layer.hpp"
#include "support/test.hpp"

namespace {

struct ColumnStorage {
  std::vector<mps::Real> theta;
  std::vector<mps::Vec3> velocity;
  std::vector<mps::Real> tracer;
  std::vector<mps::Real> mass;
  std::vector<mps::Real> exner_full;
  std::vector<mps::Real> exner_half;
  std::vector<mps::Real> height;
  std::vector<mps::Real> density_half;
  std::vector<mps::Real> momentum_k;
  std::vector<mps::Real> heat_k;
  std::vector<mps::Real> tracer_k;
  mps::Real surface_temperature = 300.0;
  mps::Real surface_exner = 1.0;
  mps::Real surface_capacity = 4.0;
  mps::Real surface_heat_conductance = 0.0;
  mps::Real surface_drag_conductance = 0.0;
  mps::Real cp = 1.0;

  [[nodiscard]] mps::BoundaryLayerColumnInput input(const mps::Real dt) const {
    return {.potential_temperature_k = theta,
            .velocity_m_s = velocity,
            .tracer_mixing_ratio = tracer,
            .air_mass_kg_m2 = mass,
            .exner_full = exner_full,
            .exner_half = exner_half,
            .height_full_m = height,
            .density_half_kg_m3 = density_half,
            .eddy_diffusivity_momentum_m2_s = momentum_k,
            .eddy_diffusivity_heat_m2_s = heat_k,
            .eddy_diffusivity_tracer_m2_s = tracer_k,
            .surface_temperature_k = surface_temperature,
            .surface_exner = surface_exner,
            .surface_heat_capacity_j_m2_k = surface_capacity,
            .surface_heat_conductance_w_m2_k = surface_heat_conductance,
            .surface_drag_conductance_kg_m2_s = surface_drag_conductance,
            .heat_capacity_cp_j_kg_k = cp,
            .time_step_s = dt};
  }
};

[[nodiscard]] ColumnStorage two_layer_column() {
  return {.theta = {2.0, 1.0},
          .velocity = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}},
          .tracer = {1.0, 0.0},
          .mass = {1.0, 1.0},
          .exner_full = {1.0, 1.0},
          .exner_half = {1.0, 1.0, 1.0},
          .height = {1.5, 0.5},
          .density_half = {1.0, 1.0, 1.0},
          .momentum_k = {0.0, 1.0, 0.0},
          .heat_k = {0.0, 1.0, 0.0},
          .tracer_k = {0.0, 1.0, 0.0}};
}

[[nodiscard]] std::pair<mps::Real, mps::Real> relaxation_step(
    const mps::Real theta, const mps::Real surface_temperature,
    const mps::Real time_step) {
  ColumnStorage column{.theta = {theta},
                       .velocity = {{0.0, 0.0, 0.0}},
                       .tracer = {0.5},
                       .mass = {2.0},
                       .exner_full = {1.0},
                       .exner_half = {1.0, 1.0},
                       .height = {0.5},
                       .density_half = {1.0, 1.0},
                       .momentum_k = {0.0, 0.0},
                       .heat_k = {0.0, 0.0},
                       .tracer_k = {0.0, 0.0},
                       .surface_temperature = surface_temperature,
                       .surface_capacity = 4.0,
                       .surface_heat_conductance = 2.0};
  const auto result = mps::implicit_boundary_layer_column(column.input(time_step));
  return {result.potential_temperature_k.front(), result.surface_temperature_k};
}

}  // namespace

MPS_TEST_CASE("fixed K backward Euler diffuses heat and tracer conservatively") {
  const auto column = two_layer_column();
  const auto result = mps::implicit_boundary_layer_column(column.input(0.25));
  MPS_CHECK_NEAR(result.potential_temperature_k[0], 11.0 / 6.0, 1e-13);
  MPS_CHECK_NEAR(result.potential_temperature_k[1], 7.0 / 6.0, 1e-13);
  MPS_CHECK_NEAR(result.tracer_mixing_ratio[0], 5.0 / 6.0, 1e-13);
  MPS_CHECK_NEAR(result.tracer_mixing_ratio[1], 1.0 / 6.0, 1e-13);
  MPS_CHECK_NEAR(result.diagnostics.heat_budget_residual_j_m2, 0.0, 1e-13);
  MPS_CHECK_NEAR(result.diagnostics.tracer_mass_change_kg_m2, 0.0, 1e-13);
  MPS_CHECK_NEAR(result.heat_flux_w_m2[1], -2.0 / 3.0, 1e-13);
}

MPS_TEST_CASE("surface and atmosphere use one implicit sensible heat flux") {
  const auto [air, surface] = relaxation_step(280.0, 300.0, 0.5);
  MPS_CHECK_NEAR(air, 2000.0 / 7.0, 1e-12);
  MPS_CHECK_NEAR(surface, 2080.0 / 7.0, 1e-12);
  ColumnStorage column{.theta = {280.0},
                       .velocity = {{0.0, 0.0, 0.0}},
                       .tracer = {0.5},
                       .mass = {2.0},
                       .exner_full = {1.0},
                       .exner_half = {1.0, 1.0},
                       .height = {0.5},
                       .density_half = {1.0, 1.0},
                       .momentum_k = {0.0, 0.0},
                       .heat_k = {0.0, 0.0},
                       .tracer_k = {0.0, 0.0},
                       .surface_temperature = 300.0,
                       .surface_capacity = 4.0,
                       .surface_heat_conductance = 2.0};
  const auto result = mps::implicit_boundary_layer_column(column.input(0.5));
  MPS_CHECK_NEAR(result.diagnostics.atmospheric_heat_change_j_m2,
                 -result.diagnostics.surface_heat_change_j_m2, 1e-12);
  MPS_CHECK_NEAR(result.diagnostics.heat_budget_residual_j_m2, 0.0, 1e-12);
}

MPS_TEST_CASE("implicit stress loss is returned as physical and numerical heat") {
  auto column = two_layer_column();
  column.theta = {300.0, 300.0};
  column.velocity = {{1.0, 0.0, 0.0}, {-1.0, 0.0, 0.0}};
  column.surface_drag_conductance = 0.5;
  const auto result = mps::implicit_boundary_layer_column(column.input(0.25));
  MPS_CHECK(result.diagnostics.physical_shear_dissipation_j_m2 > 0.0);
  MPS_CHECK(result.diagnostics.physical_surface_drag_dissipation_j_m2 > 0.0);
  MPS_CHECK(result.diagnostics.backward_euler_dissipation_j_m2 > 0.0);
  MPS_CHECK_NEAR(result.diagnostics.kinetic_energy_identity_residual_j_m2, 0.0, 1e-13);
  MPS_CHECK_NEAR(result.diagnostics.heat_budget_residual_j_m2, 0.0, 1e-12);
  MPS_CHECK_NEAR(result.diagnostics.momentum_budget_residual_kg_m_s.x, 0.0, 1e-13);
  MPS_CHECK_NEAR(result.diagnostics.momentum_budget_residual_kg_m_s.y, 0.0, 1e-13);
  MPS_CHECK_NEAR(result.diagnostics.momentum_budget_residual_kg_m_s.z, 0.0, 1e-13);
}

MPS_TEST_CASE("tracer solve obeys the maximum principle for a very large step") {
  auto column = two_layer_column();
  column.theta = {300.0, 300.0};
  column.tracer = {0.0, 1.0};
  const auto result = mps::implicit_boundary_layer_column(column.input(1e6));
  for (const auto value : result.tracer_mixing_ratio)
    MPS_CHECK(value >= 0.0 && value <= 1.0);
  MPS_CHECK_NEAR(result.diagnostics.tracer_mass_change_kg_m2, 0.0, 1e-12);
}

MPS_TEST_CASE("backward Euler surface relaxation is first order in time") {
  constexpr mps::Real exact_mean = (2.0 * 280.0 + 4.0 * 300.0) / 6.0;
  constexpr mps::Real initial_difference = 20.0;
  constexpr mps::Real decay_rate = 2.0 * (1.0 / 2.0 + 1.0 / 4.0);
  const mps::Real exact_air =
      exact_mean - (4.0 / 6.0) * initial_difference * std::exp(-decay_rate);
  std::array<mps::Real, 3> errors{};
  for (std::size_t refinement = 0; refinement < errors.size(); ++refinement) {
    const std::size_t steps = 5U << refinement;
    const mps::Real dt = 1.0 / static_cast<mps::Real>(steps);
    mps::Real air = 280.0;
    mps::Real surface = 300.0;
    for (std::size_t step = 0; step < steps; ++step) {
      const auto next = relaxation_step(air, surface, dt);
      air = next.first;
      surface = next.second;
    }
    errors[refinement] = std::abs(air - exact_air);
  }
  MPS_CHECK(std::log2(errors[0] / errors[1]) >= 0.8);
  MPS_CHECK(std::log2(errors[1] / errors[2]) >= 0.8);
}

MPS_TEST_CASE("fixed K cosine diffusion is second order in space") {
  std::array<mps::Real, 3> errors{};
  constexpr mps::Real final_time = 1e-3;
  constexpr mps::Real dt = 1e-6;
  const std::size_t steps = static_cast<std::size_t>(final_time / dt);
  for (std::size_t refinement = 0; refinement < errors.size(); ++refinement) {
    const std::size_t levels = 10U << refinement;
    const mps::Real spacing = 1.0 / static_cast<mps::Real>(levels);
    ColumnStorage column;
    column.theta.assign(levels, 300.0);
    column.velocity.assign(levels, {});
    column.tracer.resize(levels);
    column.mass.assign(levels, spacing);
    column.exner_full.assign(levels, 1.0);
    column.exner_half.assign(levels + 1, 1.0);
    column.height.resize(levels);
    column.density_half.assign(levels + 1, 1.0);
    column.momentum_k.assign(levels + 1, 0.0);
    column.heat_k.assign(levels + 1, 0.0);
    column.tracer_k.assign(levels + 1, 1.0);
    column.tracer_k.front() = 0.0;
    column.tracer_k.back() = 0.0;
    for (std::size_t level = 0; level < levels; ++level) {
      const mps::Real x = (static_cast<mps::Real>(level) + 0.5) * spacing;
      column.height[level] = 1.0 - x;
      column.tracer[level] = std::cos(std::numbers::pi_v<mps::Real> * x);
    }
    mps::BoundaryLayerColumnResult result;
    mps::BoundaryLayerColumnWorkspace workspace;
    for (std::size_t step = 0; step < steps; ++step) {
      mps::implicit_boundary_layer_column(column.input(dt), result, workspace);
      column.tracer = result.tracer_mixing_ratio;
    }
    mps::Real squared_error = 0.0;
    for (std::size_t level = 0; level < levels; ++level) {
      const mps::Real x = (static_cast<mps::Real>(level) + 0.5) * spacing;
      const mps::Real exact = std::exp(-std::numbers::pi_v<mps::Real> *
                                       std::numbers::pi_v<mps::Real> * final_time) *
                              std::cos(std::numbers::pi_v<mps::Real> * x);
      const mps::Real difference = column.tracer[level] - exact;
      squared_error += spacing * difference * difference;
    }
    errors[refinement] = std::sqrt(squared_error);
  }
  MPS_CHECK(std::log2(errors[0] / errors[1]) >= 1.7);
  MPS_CHECK(std::log2(errors[1] / errors[2]) >= 1.7);
}

int main() { return mps::test::run_all(); }
