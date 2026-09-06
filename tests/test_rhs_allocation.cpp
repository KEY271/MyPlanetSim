#include <atomic>
#include <cstdlib>
#include <new>

#include "myplanetsim/dynamics/dry_hydrostatic_driver.hpp"
#include "myplanetsim/dynamics/dry_hydrostatic_semi_implicit.hpp"
#include "support/test.hpp"

namespace {

std::atomic<bool> count_allocations{false};
std::atomic<std::size_t> allocation_count{0};

[[nodiscard]] mps::ExperimentConfig config() {
  return {.kind = mps::ExperimentKind::kDryHydrostatic,
          .planet = {6371220, 7.29212e-5, 9.80616, 287, 1004, 100000},
          .run = {0, 10, 10, 0},
          .grid = {4},
          .vertical = {.levels = 4,
                       .a_half_pa = {1000, 750, 500, 250, 0},
                       .b_half = {0, .25, .5, .75, 1},
                       .surface_pressure_pa = 100000,
                       .minimum_surface_pressure_pa = 90000,
                       .maximum_surface_pressure_pa = 110000,
                       .minimum_pressure_thickness_pa = 100,
                       .initial_temperature_k = 280,
                       .initial_potential_temperature_k = 300,
                       .temperature_floor_k = 100,
                       .transport_scheme = mps::VerticalTransportScheme::kLinear,
                       .limiter = mps::VerticalLimiterKind::kMinmod,
                       .cfl = .5},
          .dry_hydrostatic = {.test_case = mps::DryHydrostaticTestCase::kHeldSuarez},
          .physics = {.kind = mps::PhysicsKind::kHeldSuarez},
          .diagnostics = {1},
          .output_directory = "x"};
}

}  // namespace

void* operator new(const std::size_t size) {
  if (count_allocations.load(std::memory_order_relaxed))
    allocation_count.fetch_add(1, std::memory_order_relaxed);
  if (void* pointer = std::malloc(size)) return pointer;
  throw std::bad_alloc();
}

void* operator new[](const std::size_t size) { return ::operator new(size); }

void operator delete(void* pointer) noexcept { std::free(pointer); }
void operator delete[](void* pointer) noexcept { std::free(pointer); }
void operator delete(void* pointer, std::size_t) noexcept { std::free(pointer); }
void operator delete[](void* pointer, std::size_t) noexcept { std::free(pointer); }

MPS_TEST_CASE("steady dry RHS performs no heap allocations") {
  mps::DryHydrostaticDriver driver(config());
  const auto state = driver.initial_state();
  mps::DryHydrostaticRhs rhs;
  driver.rhs(state, rhs);
  driver.rhs(state, rhs);

  allocation_count.store(0, std::memory_order_relaxed);
  count_allocations.store(true, std::memory_order_relaxed);
  driver.rhs(state, rhs);
  count_allocations.store(false, std::memory_order_relaxed);

  MPS_CHECK_EQ(allocation_count.load(std::memory_order_relaxed), 0U);
}

MPS_TEST_CASE("warmed modal solve performs no heap allocations") {
  constexpr mps::PlanetParameters planet{6371220.0, 7.29212e-5, 9.80616,
                                         287.0,     1004.0,     100000.0};
  const mps::AtmosphericHybridCoordinate coordinate(
      {.a_half_pa = {1000.0, 750.0, 500.0, 250.0, 0.0},
       .b_half = {0.0, 0.25, 0.5, 0.75, 1.0}},
      80000.0, 120000.0, 100.0);
  constexpr mps::SemiImplicitParameters parameters{
      .reference_surface_pressure_pa = 100000.0,
      .reference_temperature_k = 280.0,
      .implicit_weight = 0.5,
      .wave_cfl_threshold = 0.45,
      .maximum_implicit_modes = 4,
      .nonlinear_iterations = 2,
      .linear_relative_tolerance = 1.0e-8,
      .linear_absolute_tolerance = 1.0e-12,
      .linear_maximum_iterations = 40,
      .gmres_restart = 20,
      .minimum_time_step_s = 60.0};
  const auto reference =
      mps::make_dry_hydrostatic_reference_column(coordinate, planet, parameters);
  const auto fast =
      mps::make_dry_hydrostatic_fast_operator(coordinate, planet, reference);
  const auto modes = mps::make_dry_hydrostatic_vertical_modes(reference, planet, fast);
  const mps::CubedSphereGrid grid(3, planet.radius_m);
  const auto cells = grid.cell_count();
  const auto volume = cells * fast.levels;
  mps::DryHydrostaticFastPerturbation right_hand_side{
      .surface_pressure_pa = std::vector<mps::Real>(cells),
      .horizontal_momentum_mass_kg_m_s = std::vector<mps::Vec3>(volume),
      .potential_temperature_mass_k_kg_m2 = std::vector<mps::Real>(volume)};
  for (std::size_t cell = 0; cell < cells; ++cell) {
    right_hand_side.surface_pressure_pa[cell] = 2.0 * grid.cells()[cell].center.x;
    for (std::size_t level = 0; level < fast.levels; ++level) {
      const auto n = mps::dry_hydrostatic_offset(cell, level, fast.levels);
      right_hand_side.horizontal_momentum_mass_kg_m_s[n] = mps::project_tangent(
          {grid.cells()[cell].center.y, 0.2, -grid.cells()[cell].center.x},
          grid.cells()[cell].center);
    }
  }
  const std::vector<std::size_t> selected = {0, 1};
  const mps::GmresOptions options{.restart = 20,
                                  .maximum_iterations = 40,
                                  .relative_tolerance = 1.0e-8,
                                  .absolute_tolerance = 1.0e-12};
  mps::DryHydrostaticSemiImplicitWorkspace workspace;
  mps::DryHydrostaticFastPerturbation correction;
  const auto warmup = mps::solve_dry_hydrostatic_modal_correction(
      grid, planet, fast, modes, selected, 600.0, right_hand_side, options, correction,
      workspace);
  MPS_CHECK(warmup.all_converged);

  allocation_count.store(0, std::memory_order_relaxed);
  count_allocations.store(true, std::memory_order_relaxed);
  const auto measured = mps::solve_dry_hydrostatic_modal_correction(
      grid, planet, fast, modes, selected, 600.0, right_hand_side, options, correction,
      workspace);
  count_allocations.store(false, std::memory_order_relaxed);
  MPS_CHECK(measured.all_converged);
  MPS_CHECK_EQ(allocation_count.load(std::memory_order_relaxed), 0U);
}

int main() { return mps::test::run_all(); }
