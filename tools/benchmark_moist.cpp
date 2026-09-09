#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <new>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#if defined(__unix__) || defined(__APPLE__)
#include <sys/resource.h>
#endif

#include "myplanetsim/diagnostics/dry_hydrostatic_diagnostics.hpp"
#include "myplanetsim/diagnostics/moist_diagnostics.hpp"
#include "myplanetsim/dynamics/dry_hydrostatic_driver.hpp"
#include "myplanetsim/dynamics/tracer_registry.hpp"
#include "myplanetsim/io/checkpoint.hpp"
#include "myplanetsim/io/run_metadata.hpp"
#include "myplanetsim/physics/moist_thermodynamics.hpp"
#include "myplanetsim/vertical/hybrid_pressure_coordinate.hpp"

namespace {

std::atomic<bool> count_allocations{false};
std::atomic<bool> count_warm_allocations{false};
std::atomic<std::size_t> allocation_count{0};
std::atomic<std::size_t> warm_allocation_count{0};

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

[[nodiscard]] std::vector<double> tracer_masses(const mps::CubedSphereGrid& grid,
                                                const mps::DryHydrostaticState& state,
                                                const std::size_t levels) {
  std::vector<double> result(state.tracer_count, 0.0);
  for (std::size_t tracer = 0; tracer < state.tracer_count; ++tracer)
    for (std::size_t cell = 0; cell < grid.cell_count(); ++cell)
      for (std::size_t level = 0; level < levels; ++level)
        result[tracer] += grid.cells()[cell].area_m2 *
                          state.tracer_mass_kg_m2[mps::dry_hydrostatic_tracer_offset(
                              tracer, cell, level, grid.cell_count(), levels)];
  return result;
}

[[nodiscard]] double land_water(const mps::DryHydrostaticDriver& driver,
                                const mps::DryHydrostaticState& state) {
  if (state.land_water_kg_m2.empty()) return 0.0;
  double result = 0.0;
  const auto fraction = driver.surface_boundary()->land_fraction();
  for (std::size_t cell = 0; cell < driver.grid().cell_count(); ++cell)
    result += driver.grid().cells()[cell].area_m2 * fraction[cell] *
              state.land_water_kg_m2[cell];
  return result;
}

struct StateSummary {
  double atmospheric_water_kg = 0.0;
  double land_water_kg = 0.0;
  double mean_temperature_k = 0.0;
  double mean_surface_temperature_k = 0.0;
  double maximum_wind_m_s = 0.0;
  double minimum_vapor_mixing_ratio = 0.0;
  double maximum_vapor_mixing_ratio = 0.0;
  double maximum_relative_humidity = 0.0;
};

[[nodiscard]] StateSummary summarize_state(
    const mps::DryHydrostaticDriver& driver, const mps::DryHydrostaticState& state,
    const mps::DryHydrostaticDerived& derived,
    const std::optional<std::size_t> water_vapor_tracer,
    const mps::DiluteMoistThermodynamics& thermodynamics) {
  StateSummary result;
  const std::size_t cells = driver.grid().cell_count();
  const std::size_t levels = derived.levels;
  double air_mass = 0.0;
  result.land_water_kg = land_water(driver, state);
  result.minimum_vapor_mixing_ratio =
      water_vapor_tracer ? std::numeric_limits<double>::infinity() : 0.0;
  for (std::size_t cell = 0; cell < cells; ++cell) {
    const double area = driver.grid().cells()[cell].area_m2;
    result.mean_surface_temperature_k += area * state.surface_temperature_k[cell];
    for (std::size_t level = 0; level < levels; ++level) {
      const std::size_t scalar = mps::dry_hydrostatic_offset(cell, level, levels);
      const double mass = area * derived.air_mass_kg_m2[scalar];
      air_mass += mass;
      result.mean_temperature_k += mass * derived.temperature_k[scalar];
      result.maximum_wind_m_s =
          std::max(result.maximum_wind_m_s, mps::norm(derived.velocity_m_s[scalar]));
      if (water_vapor_tracer) {
        const std::size_t vapor = mps::dry_hydrostatic_tracer_offset(
            *water_vapor_tracer, cell, level, cells, levels);
        const double mixing_ratio = derived.tracer_mixing_ratio[vapor];
        result.atmospheric_water_kg += area * state.tracer_mass_kg_m2[vapor];
        result.minimum_vapor_mixing_ratio =
            std::min(result.minimum_vapor_mixing_ratio, mixing_ratio);
        result.maximum_vapor_mixing_ratio =
            std::max(result.maximum_vapor_mixing_ratio, mixing_ratio);
        result.maximum_relative_humidity = std::max(
            result.maximum_relative_humidity,
            mps::relative_humidity(mixing_ratio, derived.temperature_k[scalar],
                                   derived.pressure_pa[scalar], thermodynamics));
      }
    }
  }
  result.mean_temperature_k /= air_mass;
  result.mean_surface_temperature_k /= driver.grid().total_area_m2();
  return result;
}

}  // namespace

void* operator new(const std::size_t size) {
  if (count_allocations.load(std::memory_order_relaxed))
    allocation_count.fetch_add(1, std::memory_order_relaxed);
  if (count_warm_allocations.load(std::memory_order_relaxed))
    warm_allocation_count.fetch_add(1, std::memory_order_relaxed);
  if (void* pointer = std::malloc(size)) return pointer;
  throw std::bad_alloc();
}

void* operator new[](const std::size_t size) { return ::operator new(size); }
void operator delete(void* pointer) noexcept { std::free(pointer); }
void operator delete[](void* pointer) noexcept { std::free(pointer); }
void operator delete(void* pointer, std::size_t) noexcept { std::free(pointer); }
void operator delete[](void* pointer, std::size_t) noexcept { std::free(pointer); }

int main(int argc, char** argv) {
  if (argc != 10 && argc != 12) {
    std::cerr << "usage: benchmark_moist CONFIG DAYS N K DT ITERATIONS TOP_PA "
                 "REFINEMENT SNAPSHOT [INITIAL_CHECKPOINT SOURCE_CONFIG]\n"
                 "ITERATIONS=-1 selects ARK2 IMEX comparison, 0 selects SSP-RK3, "
                 "and positive values select centered semi-implicit.\n";
    return 2;
  }
  try {
    auto config = mps::load_experiment_config(argv[1]);
    config.run.end_time_s = config.run.start_time_s + std::stod(argv[2]) * 86400.0;
    config.grid.cells_per_panel = std::stoi(argv[3]);
    const int levels_value = std::stoi(argv[4]);
    if (levels_value < 2 || levels_value > 1000)
      throw std::invalid_argument("invalid level count");
    const std::size_t levels = static_cast<std::size_t>(levels_value);
    config.run.time_step_s = std::stod(argv[5]);
    const int iterations = std::stoi(argv[6]);
    if (iterations < -1) throw std::invalid_argument("invalid iteration selector");
    const double top = std::stod(argv[7]);
    const double refinement = std::stod(argv[8]);
    const auto coefficients =
        mps::surface_refined_sigma_coefficients(top, levels_value, refinement);
    config.vertical.levels = levels_value;
    config.vertical.a_half_pa = coefficients.a_half_pa;
    config.vertical.b_half = coefficients.b_half;
    if (argc == 12) {
      // Keep the fixture's exact geometry, including decimal serialization.
      // The positional top/refinement arguments must still describe that grid.
      const auto source = mps::load_experiment_config(argv[11]);
      const auto matches = [](const auto& requested, const auto& stored) {
        if (requested.size() != stored.size()) return false;
        for (std::size_t i = 0; i < requested.size(); ++i)
          if (std::abs(requested[i] - stored[i]) >
              1e-9 * std::max(1.0, std::abs(requested[i])))
            return false;
        return true;
      };
      if (!matches(coefficients.a_half_pa, source.vertical.a_half_pa) ||
          !matches(coefficients.b_half, source.vertical.b_half))
        throw std::invalid_argument("benchmark checkpoint top/refinement differ");
      config.vertical.a_half_pa = source.vertical.a_half_pa;
      config.vertical.b_half = source.vertical.b_half;
    }
    config.diagnostics.interval_steps = static_cast<std::uint64_t>(
        std::max(1.0, std::ceil(86400.0 / config.run.time_step_s)));
    if (iterations != 0) {
      const auto existing =
          config.semi_implicit.value_or(mps::SemiImplicitParameters{});
      config.dry_hydrostatic.time_integrator =
          iterations == -1 ? mps::DryHydrostaticTimeIntegrator::kArk2ImexComparison
                           : mps::DryHydrostaticTimeIntegrator::kSemiImplicit;
      config.dry_hydrostatic.advective_cfl = 0.45;
      config.semi_implicit = mps::SemiImplicitParameters{
          .reference_surface_pressure_pa = config.vertical.surface_pressure_pa,
          .reference_temperature_k = config.vertical.initial_temperature_k,
          .reference_update = mps::SemiImplicitReferenceUpdate::kFixed,
          .implicit_weight = 0.5,
          .wave_cfl_threshold = 0.45,
          .maximum_implicit_modes = std::min(levels_value, 5),
          .nonlinear_iterations = std::max(iterations, 2),
          .linear_relative_tolerance = existing.linear_relative_tolerance,
          .linear_absolute_tolerance = existing.linear_absolute_tolerance,
          .linear_maximum_iterations = existing.linear_maximum_iterations,
          .gmres_restart = existing.gmres_restart,
          .minimum_time_step_s =
              std::min(existing.minimum_time_step_s, config.run.time_step_s)};
    } else {
      config.dry_hydrostatic.time_integrator =
          mps::DryHydrostaticTimeIntegrator::kExplicitSspRk3;
      config.semi_implicit.reset();
    }
    config.validate();

    const mps::DryHydrostaticDriver driver(config);
    driver.enable_rhs_profiling(std::getenv("MPS_PROFILE_RHS") != nullptr);
    auto state = driver.initial_state();
    const mps::TracerRegistry registry = config.tracers.empty()
                                             ? mps::TracerRegistry::legacy()
                                             : mps::TracerRegistry(config.tracers);
    const auto water_vapor_tracer = registry.water_vapor_index();
    const mps::DiluteMoistThermodynamics thermodynamics{
        .gas_constant_dry_air_j_kg_k = config.planet.gas_constant_j_kg_k,
        .heat_capacity_cp_j_kg_k = config.planet.heat_capacity_cp_j_kg_k};
    if (argc == 12) {
      const auto source = mps::load_experiment_config(argv[11]);
      if (source.grid.cells_per_panel != config.grid.cells_per_panel ||
          source.vertical.levels != config.vertical.levels ||
          source.vertical.a_half_pa != config.vertical.a_half_pa ||
          source.vertical.b_half != config.vertical.b_half)
        throw std::invalid_argument(
            "benchmark checkpoint grid/vertical coordinates differ");
      // Permit only time integration and diagnostic/output controls to differ.
      // This is an explicit comparison fixture import, not production restart.
      auto compatible = config;
      compatible.run = source.run;
      compatible.diagnostics = source.diagnostics;
      compatible.semi_implicit = source.semi_implicit;
      compatible.dry_hydrostatic.time_integrator =
          source.dry_hydrostatic.time_integrator;
      compatible.dry_hydrostatic.advective_cfl = source.dry_hydrostatic.advective_cfl;
      if (mps::config_fingerprint(compatible) != mps::config_fingerprint(source))
        throw std::invalid_argument(
            "benchmark checkpoint physical configuration differs");
      if (!water_vapor_tracer)
        throw std::invalid_argument("checkpoint fixture requires dilute water");
      const auto cells = driver.grid().cell_count();
      const auto checkpoint = mps::read_checkpoint_file(
          argv[10], mps::config_fingerprint(source),
          mps::dry_hydrostatic_moist_checkpoint_layout(registry),
          cells + (4 + registry.size()) * cells * levels + 2 * cells + 6);
      state = mps::unflatten_dry_hydrostatic_moist_state(
          checkpoint.time_s, checkpoint.step, checkpoint.state, cells, levels,
          registry.size());
    }
    const double initial_time_s = state.time_s;
    const double end_time_s = initial_time_s + std::stod(argv[2]) * 86400.0;
    const double initial_evaporation = state.cumulative_evaporation_kg;
    const double initial_convective_rain = state.cumulative_convective_precipitation_kg;
    const double initial_grid_rain = state.cumulative_grid_scale_precipitation_kg;
    const double initial_runoff = state.cumulative_runoff_kg;
    const double initial_ocean = state.cumulative_ocean_water_change_kg;
    const double initial_outflow = state.cumulative_external_outflow_kg;
    const auto initial_derived = driver.diagnose(state);
    const auto initial_dry = mps::diagnose_dry_hydrostatic_budgets(
        driver.grid(), state, initial_derived, config.planet);
    const auto initial_summary = summarize_state(driver, state, initial_derived,
                                                 water_vapor_tracer, thermodynamics);
    const auto initial_tracer_mass = tracer_masses(driver.grid(), state, levels);
    const double initial_system_water = initial_summary.atmospheric_water_kg +
                                        initial_summary.land_water_kg + initial_ocean +
                                        initial_outflow;

    std::ofstream samples(std::string(argv[9]) + ".samples.csv");
    if (!samples) throw std::runtime_error("cannot open moist sample diagnostics");
    samples << "day,step,mean_temperature_k,mean_surface_temperature_k,"
               "maximum_wind_m_s,atmospheric_water_kg,land_water_kg,"
               "ocean_water_change_kg,external_outflow_kg,precipitable_water_kg_m2,"
               "maximum_vapor_mixing_ratio,maximum_relative_humidity,"
               "cumulative_evaporation_kg,cumulative_convective_precipitation_kg,"
               "cumulative_grid_scale_precipitation_kg,cumulative_runoff_kg\n"
            << std::setprecision(17);
    mps::FullRhsEvaluationCounts full_rhs{};
    std::size_t accepted_steps = 0;
    std::size_t retries = 0;
    std::size_t linear_iterations = 0;
    std::size_t split_fast_operator_evaluations = 0;
    std::size_t physics_substeps = 0;
    std::size_t physics_retries = 0;
    std::size_t boundary_layer_calls = 0;
    std::size_t convection_calls = 0;
    std::size_t sbm_diagnostic_calls = 0;
    std::size_t sbm_relaxation_calls = 0;
    double accepted_seconds = 0.0;
    double toa_net_upward_energy_j = 0.0;
    double maximum_water_residual_kg = 0.0;
    double maximum_moist_enthalpy_residual_j = 0.0;
    double maximum_temperature_increment_k = 0.0;
    double maximum_vapor_increment = 0.0;
    driver.reset_rhs_profile();
    const auto start = std::chrono::steady_clock::now();
    allocation_count.store(0, std::memory_order_relaxed);
    warm_allocation_count.store(0, std::memory_order_relaxed);
    const bool measure_allocations = std::getenv("MPS_COUNT_ALLOCATIONS") != nullptr;
    count_allocations.store(measure_allocations, std::memory_order_relaxed);
    driver.advance(
        state, end_time_s,
        [&](const mps::DryHydrostaticState& sampled,
            const mps::DryHydrostaticDerived* derived,
            const mps::DryHydrostaticStepDiagnostics& step) {
          if (step.accepted_time_step_s > 0.0) {
            ++accepted_steps;
            full_rhs.initial += step.full_rhs.initial;
            full_rhs.iteration += step.full_rhs.iteration;
            full_rhs.final += step.full_rhs.final;
            full_rhs.discarded += step.full_rhs.discarded;
            accepted_seconds += step.accepted_time_step_s;
            retries += step.retry_count;
            linear_iterations += step.linear_iterations_total;
            split_fast_operator_evaluations += step.split_fast_operator_evaluations;
            physics_substeps += step.physics_substep_count;
            physics_retries += step.physics_retry_count;
            boundary_layer_calls += step.boundary_layer_column_call_count;
            convection_calls += step.convection_column_call_count;
            sbm_diagnostic_calls += step.moisture.sbm_diagnostic_column_count;
            sbm_relaxation_calls += step.moisture.sbm_relaxation_column_count;
            toa_net_upward_energy_j += step.radiation_budget.toa_net_upward_energy_j;
            maximum_water_residual_kg =
                std::max(maximum_water_residual_kg,
                         std::abs(step.moisture.water_budget_residual_kg));
            maximum_moist_enthalpy_residual_j =
                std::max(maximum_moist_enthalpy_residual_j,
                         std::abs(step.moisture.moist_enthalpy_budget_residual_j));
            maximum_temperature_increment_k =
                std::max(maximum_temperature_increment_k,
                         step.moisture.maximum_temperature_increment_k);
            maximum_vapor_increment = std::max(maximum_vapor_increment,
                                               step.moisture.maximum_vapor_increment);
          }
          if (derived != nullptr) {
            const auto summary = summarize_state(driver, sampled, *derived,
                                                 water_vapor_tracer, thermodynamics);
            samples << sampled.time_s / 86400.0 << ',' << sampled.step << ','
                    << summary.mean_temperature_k << ','
                    << summary.mean_surface_temperature_k << ','
                    << summary.maximum_wind_m_s << ',' << summary.atmospheric_water_kg
                    << ',' << summary.land_water_kg << ','
                    << sampled.cumulative_ocean_water_change_kg << ','
                    << sampled.cumulative_external_outflow_kg << ','
                    << summary.atmospheric_water_kg / driver.grid().total_area_m2()
                    << ',' << summary.maximum_vapor_mixing_ratio << ','
                    << summary.maximum_relative_humidity << ','
                    << sampled.cumulative_evaporation_kg << ','
                    << sampled.cumulative_convective_precipitation_kg << ','
                    << sampled.cumulative_grid_scale_precipitation_kg << ','
                    << sampled.cumulative_runoff_kg << '\n';
          }
          if (measure_allocations && accepted_steps == 1)
            count_warm_allocations.store(true, std::memory_order_relaxed);
        });
    count_allocations.store(false, std::memory_order_relaxed);
    count_warm_allocations.store(false, std::memory_order_relaxed);
    const double elapsed =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    if (!samples) throw std::runtime_error("cannot write moist sample diagnostics");
    if (accepted_steps == 0 || !(accepted_seconds > 0.0))
      throw std::runtime_error("moist benchmark accepted no steps");

    const auto derived = driver.diagnose(state);
    const auto final_dry = mps::diagnose_dry_hydrostatic_budgets(
        driver.grid(), state, derived, config.planet);
    const auto final_summary =
        summarize_state(driver, state, derived, water_vapor_tracer, thermodynamics);
    const auto final_tracer_mass = tracer_masses(driver.grid(), state, levels);
    const double final_system_water =
        final_summary.atmospheric_water_kg + final_summary.land_water_kg +
        state.cumulative_ocean_water_change_kg + state.cumulative_external_outflow_kg;
    const double evaporation = state.cumulative_evaporation_kg - initial_evaporation;
    const double convective_rain =
        state.cumulative_convective_precipitation_kg - initial_convective_rain;
    const double grid_rain =
        state.cumulative_grid_scale_precipitation_kg - initial_grid_rain;
    const double runoff = state.cumulative_runoff_kg - initial_runoff;
    const double exchange_scale =
        std::max({1.0, initial_system_water,
                  std::abs(evaporation) + convective_rain + grid_rain + runoff});
    double maximum_tracer_relative_drift = 0.0;
    for (std::size_t tracer = 0; tracer < final_tracer_mass.size(); ++tracer) {
      if (water_vapor_tracer && tracer == *water_vapor_tracer) continue;
      maximum_tracer_relative_drift =
          std::max(maximum_tracer_relative_drift,
                   std::abs(final_tracer_mass[tracer] - initial_tracer_mass[tracer]) /
                       std::max(1.0, std::abs(initial_tracer_mass[tracer])));
    }

    std::ofstream snapshot(argv[9]);
    if (!snapshot) throw std::runtime_error("cannot open moist snapshot");
    snapshot << "cell,level,area_m2,air_mass_kg_m2,pressure_pa,temperature_k,"
                "potential_temperature_k,wind_m_s,vapor_mixing_ratio,"
                "surface_temperature_k,land_fraction,land_water_kg_m2\n"
             << std::setprecision(17);
    const auto land_fraction = driver.surface_boundary()->land_fraction();
    for (std::size_t cell = 0; cell < driver.grid().cell_count(); ++cell)
      for (std::size_t level = 0; level < levels; ++level) {
        const std::size_t scalar = mps::dry_hydrostatic_offset(cell, level, levels);
        const double vapor =
            water_vapor_tracer
                ? derived.tracer_mixing_ratio[mps::dry_hydrostatic_tracer_offset(
                      *water_vapor_tracer, cell, level, driver.grid().cell_count(),
                      levels)]
                : 0.0;
        snapshot << cell << ',' << level << ',' << driver.grid().cells()[cell].area_m2
                 << ',' << derived.air_mass_kg_m2[scalar] << ','
                 << derived.pressure_pa[scalar] << ',' << derived.temperature_k[scalar]
                 << ',' << derived.potential_temperature_k[scalar] << ','
                 << mps::norm(derived.velocity_m_s[scalar]) << ',' << vapor << ','
                 << state.surface_temperature_k[cell] << ',' << land_fraction[cell]
                 << ','
                 << (state.land_water_kg_m2.empty() ? 0.0
                                                    : state.land_water_kg_m2[cell])
                 << '\n';
      }
    if (!snapshot) throw std::runtime_error("cannot write moist snapshot");

    const double area = driver.grid().total_area_m2();
    const double total_precipitation = convective_rain + grid_rain;
    std::cout
        << std::setprecision(17) << "elapsed_s=" << elapsed << '\n'
        << "seconds_per_model_day=" << elapsed * 86400.0 / accepted_seconds << '\n'
        << "initial_time_s=" << initial_time_s << '\n'
        << "full_rhs_initial=" << full_rhs.initial << '\n'
        << "full_rhs_iteration=" << full_rhs.iteration << '\n'
        << "full_rhs_final=" << full_rhs.final << '\n'
        << "full_rhs_discarded=" << full_rhs.discarded << '\n'
        << "full_rhs_per_step="
        << static_cast<double>(full_rhs.total()) / accepted_steps << '\n'
        << "full_rhs_per_model_day=" << full_rhs.total() * 86400.0 / accepted_seconds
        << '\n'
        << "peak_rss_bytes=" << peak_rss_bytes() << '\n'
        << "allocation_counting_enabled=" << measure_allocations << '\n'
        << "allocations=" << allocation_count.load(std::memory_order_relaxed) << '\n'
        << "allocations_per_step="
        << static_cast<double>(allocation_count.load(std::memory_order_relaxed)) /
               static_cast<double>(accepted_steps)
        << '\n'
        << "warm_allocations=" << warm_allocation_count.load(std::memory_order_relaxed)
        << '\n'
        << "warm_allocations_per_step="
        << (accepted_steps > 1 ? static_cast<double>(warm_allocation_count.load(
                                     std::memory_order_relaxed)) /
                                     static_cast<double>(accepted_steps - 1)
                               : 0.0)
        << '\n'
        << "accepted_steps=" << accepted_steps << '\n'
        << "accepted_seconds=" << accepted_seconds << '\n'
        << "retries=" << retries << '\n'
        << "linear_iterations=" << linear_iterations << '\n'
        << "split_fast_operator_evaluations=" << split_fast_operator_evaluations << '\n'
        << "physics_substeps=" << physics_substeps << '\n'
        << "physics_retries=" << physics_retries << '\n'
        << "boundary_layer_column_calls=" << boundary_layer_calls << '\n'
        << "convection_column_calls=" << convection_calls << '\n'
        << "sbm_diagnostic_column_calls=" << sbm_diagnostic_calls << '\n'
        << "sbm_relaxation_column_calls=" << sbm_relaxation_calls << '\n'
        << "dry_mass_relative_drift="
        << (final_dry.dry_mass_kg - initial_dry.dry_mass_kg) / initial_dry.dry_mass_kg
        << '\n'
        << "maximum_passive_tracer_relative_drift=" << maximum_tracer_relative_drift
        << '\n'
        << "water_budget_residual_kg=" << final_system_water - initial_system_water
        << '\n'
        << "water_budget_residual_ratio="
        << std::abs(final_system_water - initial_system_water) / exchange_scale << '\n'
        << "maximum_step_water_residual_kg=" << maximum_water_residual_kg << '\n'
        << "maximum_step_moist_enthalpy_residual_j="
        << maximum_moist_enthalpy_residual_j << '\n'
        << "mean_temperature_k=" << final_summary.mean_temperature_k << '\n'
        << "mean_surface_temperature_k=" << final_summary.mean_surface_temperature_k
        << '\n'
        << "maximum_wind_m_s=" << final_summary.maximum_wind_m_s << '\n'
        << "precipitable_water_kg_m2=" << final_summary.atmospheric_water_kg / area
        << '\n'
        << "land_water_kg_m2_global=" << final_summary.land_water_kg / area << '\n'
        << "minimum_vapor_mixing_ratio=" << final_summary.minimum_vapor_mixing_ratio
        << '\n'
        << "maximum_vapor_mixing_ratio=" << final_summary.maximum_vapor_mixing_ratio
        << '\n'
        << "maximum_relative_humidity=" << final_summary.maximum_relative_humidity
        << '\n'
        << "evaporation_kg_m2_s=" << evaporation / (area * accepted_seconds) << '\n'
        << "precipitation_kg_m2_s=" << total_precipitation / (area * accepted_seconds)
        << '\n'
        << "convective_precipitation_fraction="
        << (total_precipitation > 0.0 ? convective_rain / total_precipitation : 0.0)
        << '\n'
        << "toa_net_upward_w_m2=" << toa_net_upward_energy_j / (area * accepted_seconds)
        << '\n'
        << "maximum_temperature_increment_k=" << maximum_temperature_increment_k << '\n'
        << "maximum_vapor_increment=" << maximum_vapor_increment << '\n';
    const char* regions[] = {"diagnose_setup",  "reconstruction", "flux_cfl",
                             "column_coupling", "source",         "diffusion",
                             "physics",         "result_copy"};
    for (std::size_t i = 0; i < 8; ++i)
      std::cout << "rhs_region_" << regions[i]
                << "_s=" << driver.rhs_profile().seconds[i] << '\n';
    std::cout << "profiled_rhs_calls=" << driver.rhs_profile().calls << '\n';
    return 0;
  } catch (const std::exception& error) {
    count_allocations.store(false, std::memory_order_relaxed);
    count_warm_allocations.store(false, std::memory_order_relaxed);
    std::cerr << error.what() << '\n';
    return 1;
  }
}
