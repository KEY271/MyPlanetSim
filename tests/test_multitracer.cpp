#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "myplanetsim/dynamics/dry_hydrostatic_diffusion.hpp"
#include "myplanetsim/dynamics/dry_hydrostatic_driver.hpp"
#include "myplanetsim/dynamics/tracer_registry.hpp"
#include "myplanetsim/io/checkpoint.hpp"
#include "myplanetsim/io/run_metadata.hpp"
#include "support/test.hpp"

namespace {

[[nodiscard]] mps::ExperimentConfig base_config(const bool semi_implicit = false) {
  mps::ExperimentConfig config{
      .kind = mps::ExperimentKind::kDryHydrostatic,
      .planet = {2.0, 0.0, 10.0, 287.0, 1004.0, 100000.0},
      .run = {0.0, 0.02, 0.01, 0},
      .grid = {2},
      .vertical = {.levels = 2,
                   .a_half_pa = {1000.0, 500.0, 0.0},
                   .b_half = {0.0, 0.5, 1.0},
                   .surface_pressure_pa = 100000.0,
                   .minimum_surface_pressure_pa = 90000.0,
                   .maximum_surface_pressure_pa = 110000.0,
                   .minimum_pressure_thickness_pa = 100.0,
                   .initial_temperature_k = 280.0,
                   .initial_potential_temperature_k = 300.0,
                   .temperature_floor_k = 100.0,
                   .transport_scheme = mps::VerticalTransportScheme::kDonorCell,
                   .limiter = mps::VerticalLimiterKind::kMinmod,
                   .cfl = 0.5},
      .dry_hydrostatic = {},
      .tracers = {{.name = "first", .initial_mixing_ratio = 0.2},
                  {.name = "second", .initial_mixing_ratio = 0.4},
                  {.name = "zero", .initial_mixing_ratio = 0.0},
                  {.name = "fourth", .initial_mixing_ratio = 0.7}},
      .diagnostics = {1},
      .output_directory = "x"};
  if (semi_implicit) {
    config.dry_hydrostatic.time_integrator =
        mps::DryHydrostaticTimeIntegrator::kSemiImplicit;
    config.dry_hydrostatic.advective_cfl = 0.45;
    config.semi_implicit = mps::SemiImplicitParameters{
        .reference_surface_pressure_pa = 100000.0,
        .reference_temperature_k = 280.0,
        .reference_update = mps::SemiImplicitReferenceUpdate::kFixed,
        .implicit_weight = 0.5,
        .wave_cfl_threshold = 0.45,
        .maximum_implicit_modes = 2,
        .nonlinear_iterations = 3,
        .linear_relative_tolerance = 1e-8,
        .linear_absolute_tolerance = 1e-12,
        .linear_maximum_iterations = 40,
        .gmres_restart = 10,
        .minimum_time_step_s = 1e-5};
  }
  config.validate();
  return config;
}

}  // namespace

MPS_TEST_CASE("tracer registry resolves names and a unique water role") {
  const mps::TracerRegistry registry(
      {{.name = "dust"},
       {.name = "vapor", .role = mps::TracerRole::kWaterVapor},
       {.name = "age", .require_nonnegative = false}});
  MPS_CHECK_EQ(registry.size(), 3U);
  MPS_CHECK_EQ(registry.index_of("vapor"), 1U);
  MPS_CHECK(registry.water_vapor_index().has_value());
  MPS_CHECK_EQ(*registry.water_vapor_index(), 1U);
  MPS_CHECK_THROWS_AS(mps::TracerRegistry({{.name = "q"}, {.name = "q"}}),
                      std::invalid_argument);
  MPS_CHECK_THROWS_AS(
      mps::TracerRegistry({{.name = "a", .role = mps::TracerRole::kWaterVapor},
                           {.name = "b", .role = mps::TracerRole::kWaterVapor}}),
      std::invalid_argument);
}

MPS_TEST_CASE("multi-tracer configuration is canonical without changing legacy text") {
  auto legacy = base_config();
  legacy.tracers.clear();
  const auto legacy_text = mps::canonical_config_text(legacy);
  MPS_CHECK(legacy_text.find("tracers.") == std::string::npos);

  const auto configured = base_config();
  const auto text = mps::canonical_config_text(configured);
  std::istringstream input(text);
  const auto restored = mps::parse_experiment_config(input);
  MPS_CHECK_EQ(mps::config_fingerprint(restored), mps::config_fingerprint(configured));
  MPS_CHECK_EQ(restored.tracers.size(), 4U);
  MPS_CHECK_EQ(restored.tracers[1].name, "second");
}

MPS_TEST_CASE("components reconstruct and transport independently") {
  const auto config = base_config();
  const mps::DryHydrostaticDriver driver(config);
  auto state = driver.initial_state();
  const auto volume = state.surface_pressure_pa.size() * 2U;
  MPS_CHECK_EQ(state.tracer_count, 4U);
  MPS_CHECK_EQ(state.tracer_mass_kg_m2.size(), 4U * volume);

  // Give the first component a bounded spatial pattern, make the second an exact
  // multiple, and retain an all-zero third component.
  state.tracer_mass_kg_m2[0] *= 0.5;
  for (std::size_t n = 0; n < volume; ++n) {
    state.tracer_mass_kg_m2[volume + n] = 2.0 * state.tracer_mass_kg_m2[n];
    state.tracer_mass_kg_m2[2 * volume + n] = 0.0;
  }
  const auto rhs = driver.rhs(state);
  for (std::size_t n = 0; n < volume; ++n) {
    MPS_CHECK_NEAR(rhs.tendency.tracer_mass[volume + n],
                   2.0 * rhs.tendency.tracer_mass[n], 1e-12);
    MPS_CHECK_EQ(rhs.tendency.tracer_mass[2 * volume + n], 0.0);
  }

  auto swapped = state;
  for (std::size_t n = 0; n < volume; ++n)
    std::swap(swapped.tracer_mass_kg_m2[n], swapped.tracer_mass_kg_m2[volume + n]);
  const auto swapped_rhs = driver.rhs(swapped);
  for (std::size_t n = 0; n < volume; ++n) {
    MPS_CHECK_NEAR(swapped_rhs.tendency.tracer_mass[n],
                   rhs.tendency.tracer_mass[volume + n], 1e-12);
    MPS_CHECK_NEAR(swapped_rhs.tendency.tracer_mass[volume + n],
                   rhs.tendency.tracer_mass[n], 1e-12);
  }
}

MPS_TEST_CASE("both integrators accept one two and four tracer states") {
  for (const bool semi_implicit : {false, true}) {
    for (const std::size_t count : {1U, 2U, 4U}) {
      auto config = base_config(semi_implicit);
      config.tracers.resize(count);
      config.validate();
      const mps::DryHydrostaticDriver driver(config);
      auto state = driver.initial_state();
      driver.advance(state, config.run.end_time_s);
      MPS_CHECK_EQ(state.tracer_count, count);
      MPS_CHECK_EQ(state.time_s, config.run.end_time_s);
      MPS_CHECK(std::ranges::all_of(state.tracer_mass_kg_m2, [](const mps::Real value) {
        return value >= 0.0 && std::isfinite(value);
      }));
    }
  }
}

MPS_TEST_CASE("sharp wet-dry transport remains conservative and nonnegative") {
  for (const bool semi_implicit : {false, true}) {
    auto config = base_config(semi_implicit);
    config.run.end_time_s = 0.001;
    config.run.time_step_s = 0.001;
    config.tracers.resize(2);
    config.validate();
    const mps::DryHydrostaticDriver driver(config);
    auto state = driver.initial_state();
    const auto cells = driver.grid().cell_count();
    const auto levels = static_cast<std::size_t>(config.vertical.levels);
    const auto derived = driver.diagnose(state);
    std::vector<double> initial_mass(2, 0.0);
    for (std::size_t cell = 0; cell < cells; ++cell) {
      for (std::size_t level = 0; level < levels; ++level) {
        const auto n = mps::dry_hydrostatic_offset(cell, level, levels);
        state.horizontal_momentum_mass_kg_m_s[n] =
            derived.air_mass_kg_m2[n] * 0.05 *
            mps::cross(mps::Vec3{0.0, 0.0, 1.0}, driver.grid().cells()[cell].center);
        for (std::size_t tracer = 0; tracer < 2; ++tracer) {
          const auto q =
              mps::dry_hydrostatic_tracer_offset(tracer, cell, level, cells, levels);
          const bool wet = (cell % 2) == tracer;
          state.tracer_mass_kg_m2[q] = derived.air_mass_kg_m2[n] * (wet ? 0.01 : 0.0);
          initial_mass[tracer] +=
              driver.grid().cells()[cell].area_m2 * state.tracer_mass_kg_m2[q];
        }
      }
    }
    driver.advance(state, config.run.end_time_s);
    for (std::size_t tracer = 0; tracer < 2; ++tracer) {
      double final_mass = 0.0;
      for (std::size_t cell = 0; cell < cells; ++cell)
        for (std::size_t level = 0; level < levels; ++level) {
          const auto q =
              mps::dry_hydrostatic_tracer_offset(tracer, cell, level, cells, levels);
          MPS_CHECK(state.tracer_mass_kg_m2[q] >= -1e-12);
          final_mass +=
              driver.grid().cells()[cell].area_m2 * state.tracer_mass_kg_m2[q];
        }
      MPS_CHECK_NEAR(final_mass, initial_mass[tracer],
                     1e-11 * std::max(1.0, initial_mass[tracer]));
    }
  }
}

MPS_TEST_CASE("shared-edge tracer diffusion conserves variable-mass columns") {
  auto config = base_config();
  config.dry_hydrostatic.diffusion_kind = mps::DiffusionKind::kLaplacian;
  config.dry_hydrostatic.diffusion_coefficient = 0.01;
  config.tracers = {{.name = "water", .role = mps::TracerRole::kWaterVapor},
                    {.name = "fixed", .horizontal_diffusion = false}};
  config.validate();
  const mps::DryHydrostaticDriver driver(config);
  auto state = driver.initial_state();
  const auto cells = driver.grid().cell_count();
  const auto levels = static_cast<std::size_t>(config.vertical.levels);
  const auto volume = cells * levels;
  for (std::size_t cell = 0; cell < cells; ++cell) {
    state.surface_pressure_pa[cell] += 1000.0 * static_cast<double>(cell % 3);
    for (std::size_t level = 0; level < levels; ++level) {
      const auto q = mps::dry_hydrostatic_tracer_offset(0, cell, level, cells, levels);
      state.tracer_mass_kg_m2[q] =
          (cell == 0 ? 1.0 : 0.0) * state.tracer_mass_kg_m2[volume + q];
    }
  }
  const auto derived = driver.diagnose(state);
  const auto tendency = mps::dry_hydrostatic_diffusion_tendency(
      driver.grid(), derived, mps::DiffusionKind::kLaplacian, 0.01, config.tracers);
  for (std::size_t tracer = 0; tracer < 2; ++tracer) {
    double total = 0.0;
    for (std::size_t cell = 0; cell < cells; ++cell)
      for (std::size_t level = 0; level < levels; ++level) {
        const auto q =
            mps::dry_hydrostatic_tracer_offset(tracer, cell, level, cells, levels);
        total += driver.grid().cells()[cell].area_m2 * tendency.tracer_mass[q];
      }
    MPS_CHECK_NEAR(total, 0.0, 1e-12);
  }
  for (std::size_t n = volume; n < 2 * volume; ++n)
    MPS_CHECK_EQ(tendency.tracer_mass[n], 0.0);

  config.dry_hydrostatic.diffusion_kind = mps::DiffusionKind::kBiharmonic;
  MPS_CHECK_THROWS_AS(config.validate(), std::invalid_argument);
}

MPS_TEST_CASE("multi-tracer checkpoint binds ordered registry metadata") {
  const auto config = base_config();
  const mps::DryHydrostaticDriver driver(config);
  const auto state = driver.initial_state();
  const mps::TracerRegistry registry(config.tracers);
  const auto values = mps::flatten_dry_hydrostatic_multitracer_state(state, 2, false);
  const auto layout =
      mps::dry_hydrostatic_multitracer_checkpoint_layout(registry, false);
  std::stringstream stream;
  mps::write_checkpoint(stream, {.time_s = state.time_s,
                                 .step = state.step,
                                 .state = values,
                                 .config_fingerprint = "fingerprint",
                                 .layout_id = layout});
  const auto text = stream.str();
  MPS_CHECK(text.find("first:passive,second:passive") != std::string::npos);
  std::istringstream input(text);
  const auto checkpoint =
      mps::read_checkpoint(input, "fingerprint", layout, values.size());
  const auto restored = mps::unflatten_dry_hydrostatic_multitracer_state(
      checkpoint.time_s, checkpoint.step, checkpoint.state, driver.grid().cell_count(),
      2, registry.size(), false);
  MPS_CHECK(restored.tracer_mass_kg_m2 == state.tracer_mass_kg_m2);

  auto reordered = config.tracers;
  std::swap(reordered[0], reordered[1]);
  const auto wrong_layout = mps::dry_hydrostatic_multitracer_checkpoint_layout(
      mps::TracerRegistry(reordered), false);
  std::istringstream wrong_order(text);
  MPS_CHECK_THROWS_AS(
      mps::read_checkpoint(wrong_order, "fingerprint", wrong_layout, values.size()),
      std::runtime_error);
  MPS_CHECK_THROWS_AS(
      mps::unflatten_dry_hydrostatic_multitracer_state(
          0.0, 0, values, driver.grid().cell_count(), 2, registry.size() + 1, false),
      std::invalid_argument);
}

MPS_TEST_CASE("run metadata exposes ordered tracer names and roles") {
  const auto config = base_config();
  std::ostringstream output;
  mps::write_run_metadata(output, mps::make_run_metadata(config), config);
  MPS_CHECK(output.str().find("metadata.tracer_count = 4") != std::string::npos);
  MPS_CHECK(output.str().find("metadata.tracer.1.name = second") != std::string::npos);
  MPS_CHECK(output.str().find("metadata.tracer.1.role = passive") != std::string::npos);
}

int main() { return mps::test::run_all(); }
