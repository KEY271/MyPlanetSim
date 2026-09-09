#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
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
constexpr double kReferenceExplicitDiffusivityM2S = 1.0e5;
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
    driver.enable_rhs_profiling(std::getenv("MPS_PROFILE_RHS") != nullptr);
    auto state = driver.initial_state();

    mps::DryHydrostaticRhs rhs;
    driver.rhs(state, rhs);
    driver.rhs(state, rhs);
    allocation_count.store(0, std::memory_order_relaxed);
    count_allocations.store(true, std::memory_order_relaxed);
    driver.rhs(state, rhs);
    count_allocations.store(false, std::memory_order_relaxed);
    const std::size_t allocations = allocation_count.load(std::memory_order_relaxed);
    const auto initial_derived = driver.diagnose(state);

    const std::uint64_t initial_step = state.step;
    const double end_time =
        state.time_s + static_cast<double>(requested_steps) * config.run.time_step_s;
    driver.reset_rhs_profile();
    const auto start = std::chrono::steady_clock::now();
    driver.advance(state, end_time);
    const auto stop = std::chrono::steady_clock::now();
    mps::DryHydrostaticRhs final_rhs;
    driver.rhs(state, final_rhs);

    const std::uint64_t completed_steps = state.step - initial_step;
    const double elapsed_s = std::chrono::duration<double>(stop - start).count();
    const double seconds_per_step = elapsed_s / static_cast<double>(completed_steps);
    const std::size_t cell_levels =
        driver.grid().cell_count() * static_cast<std::size_t>(config.vertical.levels);
    const double updates_per_second =
        static_cast<double>(cell_levels) / seconds_per_step;
    const std::size_t rss_bytes = peak_rss_bytes();
    const double gamma =
        config.planet.heat_capacity_cp_j_kg_k /
        (config.planet.heat_capacity_cp_j_kg_k - config.planet.gas_constant_j_kg_k);
    const double initial_wave_speed_m_s =
        std::sqrt(gamma * config.planet.gas_constant_j_kg_k *
                  config.vertical.initial_temperature_k);
    double fast_proxy_sum_m2_s = 0.0;
    double jump_diffusivity_sum_m2_s = 0.0;
    double jump_diffusivity_min_m2_s = std::numeric_limits<double>::infinity();
    double jump_diffusivity_max_m2_s = 0.0;
    std::size_t edge_levels = 0;
    for (const auto& edge : driver.grid().edge_cache()) {
      fast_proxy_sum_m2_s += 0.5 * initial_wave_speed_m_s * edge.center_distance_m;
      for (std::size_t level = 0; level < initial_derived.levels; ++level) {
        const auto left =
            mps::dry_hydrostatic_offset(edge.left_cell, level, initial_derived.levels);
        const auto right =
            mps::dry_hydrostatic_offset(edge.right_cell, level, initial_derived.levels);
        const double speed = std::max(
            std::abs(mps::dot(initial_derived.velocity_m_s[left], edge.normal)),
            std::abs(mps::dot(initial_derived.velocity_m_s[right], edge.normal)));
        const double diffusivity = 0.5 * speed * edge.center_distance_m;
        jump_diffusivity_sum_m2_s += diffusivity;
        jump_diffusivity_min_m2_s = std::min(jump_diffusivity_min_m2_s, diffusivity);
        jump_diffusivity_max_m2_s = std::max(jump_diffusivity_max_m2_s, diffusivity);
        ++edge_levels;
      }
    }
    const double fast_proxy_mean_m2_s =
        fast_proxy_sum_m2_s / static_cast<double>(driver.grid().edge_count());
    const double jump_diffusivity_mean_m2_s =
        jump_diffusivity_sum_m2_s / static_cast<double>(edge_levels);
    const bool target_met = updates_per_second >= kMinimumCellLevelUpdatesPerSecond &&
                            seconds_per_step <= kMaximumSecondsPerStep &&
                            allocations == 0;

    std::cout
        << std::setprecision(10) << "benchmark = phase7_held_suarez\n"
        << "steps = " << completed_steps << '\n'
        << "cell_levels = " << cell_levels << '\n'
        << "elapsed_s = " << elapsed_s << '\n'
        << "seconds_per_step = " << seconds_per_step << '\n'
        << "cell_level_updates_per_s = " << updates_per_second << '\n'
        << "peak_rss_bytes = " << rss_bytes << '\n'
        << "peak_rss_bytes_per_cell_level = "
        << static_cast<double>(rss_bytes) / static_cast<double>(cell_levels) << '\n'
        << "initial_lamb_cfl_speed_m_s = " << initial_wave_speed_m_s << '\n'
        << "initial_fast_wave_stable_dt_s = "
        << rhs.horizontal_fast_wave_stable_time_step_s << '\n'
        << "initial_advective_stable_dt_s = "
        << rhs.horizontal_advective_stable_time_step_s << '\n'
        << "initial_vertical_stable_dt_s = " << rhs.vertical_stable_time_step_s << '\n'
        << "initial_diffusion_stable_dt_s = " << rhs.diffusion_stable_time_step_s
        << '\n'
        << "initial_surface_stable_dt_s = " << rhs.surface_stable_time_step_s << '\n'
        << "final_fast_wave_stable_dt_s = "
        << final_rhs.horizontal_fast_wave_stable_time_step_s << '\n'
        << "final_advective_stable_dt_s = "
        << final_rhs.horizontal_advective_stable_time_step_s << '\n'
        << "final_vertical_stable_dt_s = " << final_rhs.vertical_stable_time_step_s
        << '\n'
        << "final_diffusion_stable_dt_s = " << final_rhs.diffusion_stable_time_step_s
        << '\n'
        << "final_surface_stable_dt_s = " << final_rhs.surface_stable_time_step_s
        << '\n'
        << "inactive_fast_speed_diffusivity_proxy_m2_s = " << fast_proxy_mean_m2_s
        << '\n'
        << "initial_jump_diffusivity_min_m2_s = " << jump_diffusivity_min_m2_s << '\n'
        << "initial_jump_diffusivity_mean_m2_s = " << jump_diffusivity_mean_m2_s << '\n'
        << "initial_jump_diffusivity_max_m2_s = " << jump_diffusivity_max_m2_s << '\n'
        << "reference_explicit_diffusivity_m2_s = " << kReferenceExplicitDiffusivityM2S
        << '\n'
        << "initial_jump_to_explicit_diffusivity_ratio = "
        << jump_diffusivity_mean_m2_s / kReferenceExplicitDiffusivityM2S << '\n'
        << "allocations_per_rhs = " << allocations << '\n'
        << "registered_min_cell_level_updates_per_s = "
        << kMinimumCellLevelUpdatesPerSecond << '\n'
        << "registered_max_seconds_per_step = " << kMaximumSecondsPerStep << '\n'
        << "registered_target_met = " << std::boolalpha << target_met << '\n';
    const char* regions[] = {"diagnose_setup",  "reconstruction", "flux_cfl",
                             "column_coupling", "source",         "diffusion",
                             "physics",         "result_copy"};
    for (std::size_t i = 0; i < 8; ++i)
      std::cout << "rhs_region_" << regions[i]
                << "_s=" << driver.rhs_profile().seconds[i] << '\n';
    std::cout << "profiled_rhs_calls=" << driver.rhs_profile().calls << '\n'
              << "prepared_reconstruction_bytes="
              << driver.rhs_profile().prepared_reconstruction_bytes << '\n'
              << "edge_flux_bytes=" << driver.rhs_profile().edge_flux_bytes << '\n'
              << "eliminated_face_state_bytes="
              << driver.rhs_profile().eliminated_face_state_bytes << '\n';
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "benchmark error: " << error.what() << '\n';
    return 1;
  }
}
