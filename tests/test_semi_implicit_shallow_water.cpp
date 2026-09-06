#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <vector>

#include "myplanetsim/numerics/helmholtz_solver.hpp"
#include "support/test.hpp"

namespace {

struct LinearState {
  std::vector<mps::Real> surface;
  std::vector<mps::Real> edge_velocity;
};

struct Fixture {
  mps::CubedSphereGrid grid{4, 1.0};
  mps::Real gravity = 1.0;
  mps::Real depth = 1.0;
  mps::GmresWorkspace workspace;
  std::vector<mps::Real> laplacian;
  std::vector<mps::Real> divergence;
  std::vector<mps::Real> right_hand_side;
  std::vector<mps::Real> next_surface;

  Fixture()
      : laplacian(grid.cell_count()),
        divergence(grid.cell_count()),
        right_hand_side(grid.cell_count()),
        next_surface(grid.cell_count()) {}

  [[nodiscard]] mps::Real minimum_spacing() const {
    mps::Real result = std::numeric_limits<mps::Real>::infinity();
    for (const auto& edge : grid.edges())
      result = std::min(result, grid.edge_cache()[edge.id].center_distance_m);
    return result;
  }

  [[nodiscard]] LinearState initial_state() const {
    LinearState state{.surface = std::vector<mps::Real>(grid.cell_count()),
                      .edge_velocity = std::vector<mps::Real>(grid.edge_count(), 0.0)};
    mps::Real mean = 0.0;
    for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
      state.surface[cell] = 1.0e-3 * grid.cells()[cell].center.x;
      mean += grid.cells()[cell].area_m2 * state.surface[cell];
    }
    mean /= grid.total_area_m2();
    for (auto& value : state.surface) value -= mean;
    return state;
  }

  void diagnose_divergence(const std::span<const mps::Real> edge_velocity) {
    std::fill(divergence.begin(), divergence.end(), 0.0);
    for (const auto& edge : grid.edges()) {
      const auto& cached = grid.edge_cache()[edge.id];
      const mps::Real flux = edge.length_m * edge_velocity[edge.id];
      divergence[cached.left_cell] += flux / grid.cells()[cached.left_cell].area_m2;
      divergence[cached.right_cell] -= flux / grid.cells()[cached.right_cell].area_m2;
    }
  }

  [[nodiscard]] mps::GmresResult step(LinearState& state, const mps::Real time_step) {
    const mps::Real coefficient = 0.25 * time_step * time_step * gravity * depth;
    const mps::FiniteVolumeHelmholtzOperator helmholtz(grid, coefficient);
    helmholtz.apply_laplacian(state.surface, laplacian);
    diagnose_divergence(state.edge_velocity);
    for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
      right_hand_side[cell] = state.surface[cell] + coefficient * laplacian[cell] -
                              time_step * depth * divergence[cell];
    }
    next_surface = state.surface;
    const auto solve =
        mps::solve_finite_volume_helmholtz(helmholtz, right_hand_side, next_surface,
                                           {.restart = 20,
                                            .maximum_iterations = 80,
                                            .relative_tolerance = 1.0e-12,
                                            .absolute_tolerance = 1.0e-14},
                                           workspace);
    if (!solve.converged()) return solve;
    for (const auto& edge : grid.edges()) {
      const auto& cached = grid.edge_cache()[edge.id];
      const mps::Real old_gradient =
          (state.surface[cached.right_cell] - state.surface[cached.left_cell]) /
          cached.center_distance_m;
      const mps::Real new_gradient =
          (next_surface[cached.right_cell] - next_surface[cached.left_cell]) /
          cached.center_distance_m;
      state.edge_velocity[edge.id] -=
          0.5 * time_step * gravity * (old_gradient + new_gradient);
    }
    state.surface.swap(next_surface);
    return solve;
  }

  [[nodiscard]] mps::Real mass(const LinearState& state) const {
    mps::Real result = 0.0;
    for (std::size_t cell = 0; cell < grid.cell_count(); ++cell)
      result += grid.cells()[cell].area_m2 * state.surface[cell];
    return result;
  }

  [[nodiscard]] mps::Real energy(const LinearState& state) const {
    mps::Real result = 0.0;
    for (std::size_t cell = 0; cell < grid.cell_count(); ++cell)
      result += 0.5 * gravity * grid.cells()[cell].area_m2 * state.surface[cell] *
                state.surface[cell];
    for (const auto& edge : grid.edges()) {
      result += 0.5 * depth * edge.length_m *
                grid.edge_cache()[edge.id].center_distance_m *
                state.edge_velocity[edge.id] * state.edge_velocity[edge.id];
    }
    return result;
  }
};

[[nodiscard]] mps::Real state_error(const LinearState& left, const LinearState& right) {
  mps::Real sum = 0.0;
  for (std::size_t index = 0; index < left.surface.size(); ++index) {
    const mps::Real difference = left.surface[index] - right.surface[index];
    sum += difference * difference;
  }
  for (std::size_t index = 0; index < left.edge_velocity.size(); ++index) {
    const mps::Real difference = left.edge_velocity[index] - right.edge_velocity[index];
    sum += difference * difference;
  }
  return std::sqrt(sum);
}

[[nodiscard]] LinearState integrate(const mps::Real time_step,
                                    const mps::Real end_time) {
  Fixture fixture;
  auto state = fixture.initial_state();
  const auto steps = static_cast<std::size_t>(std::llround(end_time / time_step));
  for (std::size_t step = 0; step < steps; ++step)
    MPS_CHECK(fixture.step(state, time_step).converged());
  return state;
}

}  // namespace

MPS_TEST_CASE("linear shallow-water CN is stable through wave Courant eight") {
  for (const mps::Real courant : std::array<mps::Real, 5>{0.5, 1.0, 2.0, 4.0, 8.0}) {
    Fixture fixture;
    auto state = fixture.initial_state();
    const mps::Real initial_mass = fixture.mass(state);
    const mps::Real initial_energy = fixture.energy(state);
    const mps::Real time_step = courant * fixture.minimum_spacing();
    std::size_t maximum_iterations = 0;
    for (int step = 0; step < 20; ++step) {
      const auto solve = fixture.step(state, time_step);
      MPS_CHECK(solve.converged());
      maximum_iterations = std::max(maximum_iterations, solve.iterations);
    }
    MPS_CHECK(maximum_iterations < 80);
    MPS_CHECK_NEAR(fixture.mass(state), initial_mass, 2.0e-12);
    MPS_CHECK_NEAR(fixture.energy(state), initial_energy, 2.0e-9 * initial_energy);
  }
}

MPS_TEST_CASE("linear shallow-water CN has second-order temporal convergence") {
  Fixture scale;
  const mps::Real spacing = scale.minimum_spacing();
  const mps::Real end_time = 6.4 * spacing;
  const auto reference = integrate(0.025 * spacing, end_time);
  const auto coarse = integrate(0.4 * spacing, end_time);
  const auto medium = integrate(0.2 * spacing, end_time);
  const auto fine = integrate(0.1 * spacing, end_time);
  const mps::Real coarse_error = state_error(coarse, reference);
  const mps::Real medium_error = state_error(medium, reference);
  const mps::Real fine_error = state_error(fine, reference);
  MPS_CHECK(coarse_error / medium_error > 3.5);
  MPS_CHECK(medium_error / fine_error > 3.5);
}

int main() { return mps::test::run_all(); }
