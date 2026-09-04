#include <atomic>
#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <new>
#include <stdexcept>
#include <string_view>

#if defined(__unix__) || defined(__APPLE__)
#include <sys/resource.h>
#endif

#include "myplanetsim/config/experiment_config.hpp"
#include "myplanetsim/dynamics/dry_hydrostatic_driver.hpp"

namespace {

constexpr double kMinimumCellLevelUpdatesPerSecond = 8.0e5;
constexpr double kMaximumSecondsPerStep = 0.022;
std::atomic<bool> count_allocations{false};
std::atomic<std::size_t> allocation_count{0};

[[nodiscard]] std::size_t parse_steps(const int argc, char** argv) {
  if (argc == 1) return 100;
  if (argc != 3 || std::string_view(argv[1]) != "--steps")
    throw std::invalid_argument("usage: benchmark_dry_core [--steps COUNT]");
  const std::string_view text(argv[2]);
  std::size_t value = 0;
  for (const char character : text) {
    if (character < '0' || character > '9')
      throw std::invalid_argument("benchmark step count must be a positive integer");
    value = value * 10 + static_cast<std::size_t>(character - '0');
  }
  if (value == 0) throw std::invalid_argument("benchmark step count must be positive");
  return value;
}

[[nodiscard]] std::size_t peak_rss_bytes() noexcept {
#if defined(__unix__) || defined(__APPLE__)
  rusage usage{};
  if (getrusage(RUSAGE_SELF, &usage) != 0) return 0;
#if defined(__APPLE__)
  return static_cast<std::size_t>(usage.ru_maxrss);
#else
  return static_cast<std::size_t>(usage.ru_maxrss) * 1024;
#endif
#else
  return 0;
#endif
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

int main(int argc, char** argv) {
  try {
    const std::size_t requested_steps = parse_steps(argc, argv);
    auto config = mps::load_experiment_config(MPS_DRY_CORE_BENCHMARK_CONFIG);
    const mps::DryHydrostaticDriver driver(config);
    auto state = driver.initial_state();

    mps::DryHydrostaticRhs rhs;
    driver.rhs(state, rhs);
    driver.rhs(state, rhs);
    allocation_count.store(0, std::memory_order_relaxed);
    count_allocations.store(true, std::memory_order_relaxed);
    driver.rhs(state, rhs);
    count_allocations.store(false, std::memory_order_relaxed);
    const std::size_t allocations = allocation_count.load(std::memory_order_relaxed);

    const std::uint64_t initial_step = state.step;
    const double end_time =
        state.time_s + static_cast<double>(requested_steps) * config.run.time_step_s;
    const auto start = std::chrono::steady_clock::now();
    driver.advance(state, end_time);
    const auto stop = std::chrono::steady_clock::now();

    const std::uint64_t completed_steps = state.step - initial_step;
    const double elapsed_s = std::chrono::duration<double>(stop - start).count();
    const double seconds_per_step = elapsed_s / static_cast<double>(completed_steps);
    const std::size_t cell_levels =
        driver.grid().cell_count() * static_cast<std::size_t>(config.vertical.levels);
    const double updates_per_second =
        static_cast<double>(cell_levels) / seconds_per_step;
    const std::size_t rss_bytes = peak_rss_bytes();
    const bool target_met = updates_per_second >= kMinimumCellLevelUpdatesPerSecond &&
                            seconds_per_step <= kMaximumSecondsPerStep &&
                            allocations == 0;

    std::cout << std::setprecision(10) << "benchmark = phase7_held_suarez\n"
              << "steps = " << completed_steps << '\n'
              << "cell_levels = " << cell_levels << '\n'
              << "elapsed_s = " << elapsed_s << '\n'
              << "seconds_per_step = " << seconds_per_step << '\n'
              << "cell_level_updates_per_s = " << updates_per_second << '\n'
              << "peak_rss_bytes = " << rss_bytes << '\n'
              << "peak_rss_bytes_per_cell_level = "
              << static_cast<double>(rss_bytes) / static_cast<double>(cell_levels)
              << '\n'
              << "allocations_per_rhs = " << allocations << '\n'
              << "registered_min_cell_level_updates_per_s = "
              << kMinimumCellLevelUpdatesPerSecond << '\n'
              << "registered_max_seconds_per_step = " << kMaximumSecondsPerStep << '\n'
              << "registered_target_met = " << std::boolalpha << target_met << '\n';
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "benchmark error: " << error.what() << '\n';
    return 1;
  }
}
