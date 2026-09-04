#include <atomic>
#include <cstdlib>
#include <new>

#include "myplanetsim/dynamics/dry_hydrostatic_driver.hpp"
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

int main() { return mps::test::run_all(); }
