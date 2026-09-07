// Phase 12 P12.07 single-column experiments. The driver couplings are gated by
// test_dry_mixing_core.cpp; here the same first-order splitting is reproduced on one
// column so heating, nocturnal cooling, the radiative-convective balance, and the
// vertical/temporal sensitivity can be measured without global dynamics.
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "myplanetsim/physics/boundary_layer.hpp"
#include "myplanetsim/physics/dry_convective_adjustment.hpp"
#include "myplanetsim/physics/gray_radiation.hpp"
#include "myplanetsim/vertical/hybrid_pressure_coordinate.hpp"
#include "myplanetsim/vertical/hydrostatic_column.hpp"
#include "support/test.hpp"

namespace {

constexpr mps::Real kGravity = 9.80616;
constexpr mps::Real kGasConstant = 287.0;
constexpr mps::Real kHeatCapacityCp = 1004.0;
constexpr mps::Real kReferencePressure = 100000.0;
constexpr mps::Real kSurfacePressure = 100000.0;
constexpr mps::Real kTopPressure = 1000.0;
constexpr mps::Real kSurfaceRefinement = 5.0;
constexpr mps::Real kDay = 86400.0;

[[nodiscard]] mps::RadiationParameters radiation_parameters() {
  return mps::RadiationParameters{
      .shortwave_absorption_m2_kg = 2e-6,
      .longwave_absorption_ref_m2_kg = 3e-4,
      .reference_pressure_pa = kReferencePressure,
      .longwave_pressure_exponent = 1.0,
      .longwave_diffusivity_factor = 1.66,
      .shortwave_diffuse_factor = 1.66,
      .cfl = 0.5,
  };
}

struct ColumnOptions {
  std::size_t levels = 20;
  mps::Real surface_temperature_k = 300.0;
  mps::Real surface_heat_capacity_j_m2_k = 1e6;
  mps::Real bottom_potential_temperature_k = 300.0;
  mps::Real potential_temperature_lapse_k_m = 4e-3;
  mps::Real wind_m_s = 10.0;
  mps::Real stellar_flux_w_m2 = 0.0;
  bool radiation = false;
  bool diurnal = false;
  mps::Real gustiness_m_s = 1.0;
  mps::Real critical_richardson = 1.0;
};

// One column advanced with the Phase 12 order: radiation, then the implicit boundary
// layer with its surface reservoir, then the dry convective adjustment.
class MixingColumn {
 public:
  explicit MixingColumn(const ColumnOptions& options)
      : options_(options),
        coordinate_(mps::surface_refined_sigma_coefficients(
                        kTopPressure, static_cast<mps::Index>(options.levels),
                        kSurfaceRefinement),
                    0.5 * kSurfacePressure, 2.0 * kSurfacePressure, 1.0) {
    surface_temperature_k = options.surface_temperature_k;
    geometry_ = coordinate_.geometry(kSurfacePressure, kGravity, kGasConstant,
                                     kHeatCapacityCp, kReferencePressure);
    potential_temperature_k.assign(options.levels,
                                   options.bottom_potential_temperature_k);
    velocity_m_s.assign(options.levels, mps::Vec3{options.wind_m_s, 0.0, 0.0});
    tracer_mixing_ratio.assign(options.levels, 1.0);
    update_heights();
    for (std::size_t level = 0; level < options.levels; ++level)
      potential_temperature_k[level] =
          options.bottom_potential_temperature_k +
          options.potential_temperature_lapse_k_m * height_full_m[level];
    update_heights();
  }

  void step(const mps::Real time_s, const mps::Real dt) {
    if (options_.radiation) apply_radiation(time_s, dt);
    update_heights();
    apply_boundary_layer(dt);
    apply_convection();
  }

  [[nodiscard]] mps::Real enthalpy_j_m2() const {
    mps::Real result = 0.0;
    for (std::size_t level = 0; level < options_.levels; ++level)
      result += kHeatCapacityCp * geometry_.air_mass_kg_m2[level] *
                geometry_.exner_full[level] * potential_temperature_k[level];
    return result;
  }

  [[nodiscard]] mps::Real kinetic_energy_j_m2() const {
    mps::Real result = 0.0;
    for (std::size_t level = 0; level < options_.levels; ++level)
      result += 0.5 * geometry_.air_mass_kg_m2[level] *
                dot(velocity_m_s[level], velocity_m_s[level]);
    return result;
  }

  [[nodiscard]] mps::Real surface_storage_j_m2() const {
    return options_.surface_heat_capacity_j_m2_k * surface_temperature_k;
  }

  // Depth of the neutral layer resting on the surface, measured from the interface
  // above the highest layer whose potential temperature still matches the bottom one.
  [[nodiscard]] mps::Real mixed_depth_m() const {
    std::size_t level = options_.levels - 1;
    while (level > 0 && std::abs(potential_temperature_k[level - 1] -
                                 potential_temperature_k[options_.levels - 1]) < 1e-6)
      --level;
    return height_half_m[level];
  }

  [[nodiscard]] const mps::HybridPressureGeometry& geometry() const {
    return geometry_;
  }
  [[nodiscard]] const std::vector<mps::Real>& height_full() const {
    return height_full_m;
  }

  std::vector<mps::Real> potential_temperature_k;
  std::vector<mps::Vec3> velocity_m_s;
  std::vector<mps::Real> tracer_mixing_ratio;
  std::vector<mps::Real> height_half_m;
  std::vector<mps::Real> height_full_m;
  mps::Real surface_temperature_k = 0.0;
  mps::BoundaryLayerBulkResult bulk;
  mps::BoundaryLayerColumnResult mixed;
  mps::DryConvectiveAdjustmentResult adjusted;
  mps::Real top_of_atmosphere_downward_w_m2 = 0.0;
  mps::Real returned_dissipation_j_m2 = 0.0;

 private:
  void update_heights() {
    const auto column = mps::integrate_hydrostatic_column(
        geometry_, potential_temperature_k, kHeatCapacityCp, kGravity, 0.0);
    height_half_m = column.height_half_m;
    height_full_m = column.height_full_m;
    height_half_m.back() = 0.0;
  }

  void apply_radiation(const mps::Real time_s, const mps::Real dt) {
    std::vector<mps::Real> temperature(options_.levels);
    for (std::size_t level = 0; level < options_.levels; ++level)
      temperature[level] = geometry_.exner_full[level] * potential_temperature_k[level];
    const mps::Real zenith =
        options_.diurnal ? std::max(0.0, std::cos(2.0 * kPi * time_s / kDay)) : 0.25;
    const mps::GrayRadiativeColumnSourceInput input{
        .radiation =
            mps::GrayRadiationColumnInput{
                .pressure_half_pa = geometry_.pressure_half_pa,
                .temperature_k = temperature,
                .surface_temperature_k = surface_temperature_k,
                .gravity_m_s2 = kGravity,
                .stellar_flux_w_m2 = options_.stellar_flux_w_m2,
                .cosine_solar_zenith = zenith,
                .surface_albedo = 0.2,
                .surface_emissivity = 1.0,
                .parameters = radiation_parameters(),
            },
        .exner_full = geometry_.exner_full,
        .heat_capacity_cp_j_kg_k = kHeatCapacityCp,
        .surface_heat_capacity_j_m2_k = options_.surface_heat_capacity_j_m2_k,
        // The boundary layer owns the sensible heat flux in Phase 12.
        .air_exchange_coefficient_w_m2_k = 0.0,
        .internal_heat_flux_w_m2 = 0.0,
    };
    const auto tendency = mps::gray_radiative_column_tendency(input);
    top_of_atmosphere_downward_w_m2 = -tendency.fluxes.net_flux_w_m2.front();
    for (std::size_t level = 0; level < options_.levels; ++level)
      potential_temperature_k[level] +=
          dt * tendency.potential_temperature_mass_k_kg_m2_s[level] /
          geometry_.air_mass_kg_m2[level];
    surface_temperature_k += dt * tendency.surface_temperature_k_s;
  }

  void apply_boundary_layer(const mps::Real dt) {
    mps::diagnose_boundary_layer_column(
        {.potential_temperature_k = potential_temperature_k,
         .temperature_k = temperature_view(),
         .velocity_m_s = velocity_m_s,
         .pressure_half_pa = geometry_.pressure_half_pa,
         .exner_half = geometry_.exner_half,
         .height_half_m = height_half_m,
         .height_full_m = height_full_m,
         .surface_temperature_k = surface_temperature_k,
         .surface_exner = geometry_.exner_half.back(),
         .gravity_m_s2 = kGravity,
         .gas_constant_j_kg_k = kGasConstant,
         .heat_capacity_cp_j_kg_k = kHeatCapacityCp,
         .critical_richardson = options_.critical_richardson,
         .turbulent_prandtl = 1.0,
         .gustiness_m_s = options_.gustiness_m_s,
         .land_fraction = 1.0,
         .land_roughness = {.momentum_m = 0.1, .heat_m = 0.01},
         .ocean_roughness = {.momentum_m = 0.1, .heat_m = 0.01}},
        bulk);
    mps::implicit_boundary_layer_column(
        {.potential_temperature_k = potential_temperature_k,
         .velocity_m_s = velocity_m_s,
         .tracer_mixing_ratio = tracer_mixing_ratio,
         .air_mass_kg_m2 = geometry_.air_mass_kg_m2,
         .exner_full = geometry_.exner_full,
         .exner_half = geometry_.exner_half,
         .height_full_m = height_full_m,
         .density_half_kg_m3 = bulk.density_half_kg_m3,
         .eddy_diffusivity_momentum_m2_s = bulk.eddy_diffusivity_momentum_m2_s,
         .eddy_diffusivity_heat_m2_s = bulk.eddy_diffusivity_heat_m2_s,
         .eddy_diffusivity_tracer_m2_s = bulk.eddy_diffusivity_tracer_m2_s,
         .surface_temperature_k = surface_temperature_k,
         .surface_exner = geometry_.exner_half.back(),
         .surface_heat_capacity_j_m2_k = options_.surface_heat_capacity_j_m2_k,
         .surface_heat_conductance_w_m2_k = bulk.surface_heat_conductance_w_m2_k,
         .surface_drag_conductance_kg_m2_s = bulk.surface_drag_conductance_kg_m2_s,
         .heat_capacity_cp_j_kg_k = kHeatCapacityCp,
         .time_step_s = dt},
        mixed, mixing_workspace_);
    potential_temperature_k = mixed.potential_temperature_k;
    velocity_m_s = mixed.velocity_m_s;
    tracer_mixing_ratio = mixed.tracer_mixing_ratio;
    surface_temperature_k = mixed.surface_temperature_k;
    returned_dissipation_j_m2 = mixed.diagnostics.returned_dissipation_heat_j_m2;
  }

  void apply_convection() {
    mps::dry_convective_adjustment(
        {.potential_temperature_k = potential_temperature_k,
         .air_mass_kg_m2 = geometry_.air_mass_kg_m2,
         .exner_full = geometry_.exner_full,
         .exner_half = geometry_.exner_half,
         .heat_capacity_cp_j_kg_k = kHeatCapacityCp,
         .heat_capacity_cv_j_kg_k = kHeatCapacityCp - kGasConstant,
         .stability_tolerance_k = 1e-10},
        adjusted, adjustment_workspace_);
    potential_temperature_k = adjusted.adjusted_potential_temperature_k;
  }

  [[nodiscard]] std::span<const mps::Real> temperature_view() {
    temperature_k_.resize(options_.levels);
    for (std::size_t level = 0; level < options_.levels; ++level)
      temperature_k_[level] =
          geometry_.exner_full[level] * potential_temperature_k[level];
    return temperature_k_;
  }

  static constexpr mps::Real kPi = 3.14159265358979323846;

  ColumnOptions options_;
  mps::AtmosphericHybridCoordinate coordinate_;
  mps::HybridPressureGeometry geometry_;
  std::vector<mps::Real> temperature_k_;
  mps::BoundaryLayerColumnWorkspace mixing_workspace_;
  mps::DryConvectiveAdjustmentWorkspace adjustment_workspace_;
};

[[nodiscard]] std::size_t layers_below(const MixingColumn& column,
                                       const mps::Real height_m) {
  std::size_t count = 0;
  for (const mps::Real height : column.height_full())
    if (height < height_m) ++count;
  return count;
}

}  // namespace

MPS_TEST_CASE("surface-refined sigma resolves the layers a boundary layer needs") {
  for (const mps::Index levels : {20, 40, 80}) {
    const auto coefficients = mps::surface_refined_sigma_coefficients(
        kTopPressure, levels, kSurfaceRefinement);
    coefficients.validate(0.5 * kSurfacePressure, 2.0 * kSurfacePressure, 1.0);
    MPS_CHECK(!mps::is_uniform_sigma(coefficients.a_half_pa, coefficients.b_half));
    MPS_CHECK_EQ(coefficients.b_half.front(), 0.0);
    MPS_CHECK_EQ(coefficients.b_half.back(), 1.0);
    for (std::size_t level = 1; level + 1 < coefficients.b_half.size(); ++level) {
      const mps::Real upper =
          coefficients.b_half[level] - coefficients.b_half[level - 1];
      const mps::Real lower =
          coefficients.b_half[level + 1] - coefficients.b_half[level];
      MPS_CHECK(lower < upper);
    }
    const MixingColumn column({.levels = static_cast<std::size_t>(levels)});
    // The registered starting target: five resolved layers in the lowest kilometre.
    MPS_CHECK(layers_below(column, 1000.0) >= 5);
  }
  // The lowest full level is reported rather than assumed; K=20 and K=40 sit inside the
  // 30--100 m target of the plan and K=80 refines below it.
  const MixingColumn coarse({.levels = 20});
  const MixingColumn medium({.levels = 40});
  const MixingColumn fine({.levels = 80});
  MPS_CHECK(coarse.height_full().back() > 30.0 && coarse.height_full().back() < 100.0);
  MPS_CHECK(medium.height_full().back() > 30.0 && medium.height_full().back() < 100.0);
  MPS_CHECK(fine.height_full().back() < medium.height_full().back());
  MPS_CHECK(fine.height_full().back() > 0.0);
}

MPS_TEST_CASE("surface refinement of one reproduces the uniform sigma ramp") {
  const auto refined = mps::surface_refined_sigma_coefficients(kTopPressure, 12, 1.0);
  const auto uniform = mps::uniform_sigma_coefficients(kTopPressure, 12);
  MPS_CHECK(refined.a_half_pa == uniform.a_half_pa);
  MPS_CHECK(refined.b_half == uniform.b_half);
  MPS_CHECK_THROWS_AS(mps::surface_refined_sigma_coefficients(kTopPressure, 12, 0.5),
                      std::invalid_argument);
}

MPS_TEST_CASE("daytime surface heating grows a neutral mixed layer") {
  MixingColumn column({.levels = 20,
                       .surface_temperature_k = 320.0,
                       .surface_heat_capacity_j_m2_k = 1e9});
  const mps::Real initial_depth = column.mixed_depth_m();
  mps::Real previous_depth = initial_depth;
  mps::Real time_s = 0.0;
  bool grew = false;
  for (int step = 0; step < 240; ++step) {
    column.step(time_s, 60.0);
    time_s += 60.0;
    MPS_CHECK(column.bulk.diagnostics.sensible_heat_flux_w_m2 > 0.0);
    MPS_CHECK(column.bulk.diagnostics.surface_richardson < 0.0);
    MPS_CHECK_EQ(column.adjusted.diagnostics.unstable_interface_fraction_after, 0.0);
    const mps::Real depth = column.mixed_depth_m();
    MPS_CHECK(depth >= previous_depth - 1e-9);
    if (depth > previous_depth) grew = true;
    previous_depth = depth;
  }
  MPS_CHECK(grew);
  MPS_CHECK(column.mixed_depth_m() > initial_depth);
  MPS_CHECK(column.potential_temperature_k.back() > 300.0);
  MPS_CHECK(column.bulk.diagnostics.boundary_layer_height_m > 0.0);
  MPS_CHECK(!column.bulk.diagnostics.reaches_model_top);
}

MPS_TEST_CASE("nocturnal cooling leaves a shallow stable layer without adjustment") {
  MixingColumn day({.levels = 20,
                    .surface_temperature_k = 320.0,
                    .surface_heat_capacity_j_m2_k = 1e9});
  MixingColumn night({.levels = 20,
                      .surface_temperature_k = 285.0,
                      .surface_heat_capacity_j_m2_k = 1e9});
  mps::Real time_s = 0.0;
  for (int step = 0; step < 120; ++step) {
    day.step(time_s, 60.0);
    night.step(time_s, 60.0);
    time_s += 60.0;
  }
  MPS_CHECK(night.bulk.diagnostics.sensible_heat_flux_w_m2 < 0.0);
  MPS_CHECK(night.bulk.diagnostics.surface_richardson > 0.0);
  MPS_CHECK(night.bulk.diagnostics.surface_stability_factor < 1.0);
  MPS_CHECK(night.bulk.diagnostics.boundary_layer_height_m <
            day.bulk.diagnostics.boundary_layer_height_m);
  // A stably cooled column has nothing for the dry adjustment to do.
  MPS_CHECK_EQ(night.adjusted.diagnostics.adjusted_layer_count, 0U);
  MPS_CHECK(night.potential_temperature_k.back() <
            night.potential_temperature_k[night.potential_temperature_k.size() - 2]);
  MPS_CHECK(night.mixed_depth_m() < day.mixed_depth_m());
}

MPS_TEST_CASE(
    "radiative-convective column closes energy against the top of atmosphere") {
  MixingColumn column({.levels = 20,
                       .surface_heat_capacity_j_m2_k = 2e6,
                       .stellar_flux_w_m2 = 1361.0,
                       .radiation = true,
                       .diurnal = true});
  const mps::Real initial_energy = column.enthalpy_j_m2() +
                                   column.kinetic_energy_j_m2() +
                                   column.surface_storage_j_m2();
  const mps::Real dt = 120.0;
  mps::Real absorbed_j_m2 = 0.0;
  mps::Real scale = 0.0;
  mps::Real time_s = 0.0;
  for (int step = 0; step < 720; ++step) {
    column.step(time_s, dt);
    time_s += dt;
    absorbed_j_m2 += dt * column.top_of_atmosphere_downward_w_m2;
    scale += dt * std::abs(column.top_of_atmosphere_downward_w_m2);
    MPS_CHECK_EQ(column.adjusted.diagnostics.unstable_interface_fraction_after, 0.0);
    MPS_CHECK(std::isfinite(column.surface_temperature_k));
  }
  const mps::Real final_energy = column.enthalpy_j_m2() + column.kinetic_energy_j_m2() +
                                 column.surface_storage_j_m2();
  // Convection conserves enthalpy, the boundary layer moves heat between the column and
  // the reservoir and returns the kinetic energy it removes, so the only net source
  // over a day is the radiative flux through the top.
  MPS_CHECK_NEAR(final_energy - initial_energy, absorbed_j_m2, 1e-9 * scale);
  MPS_CHECK(column.surface_temperature_k > 100.0);
}

MPS_TEST_CASE("coupled column splitting is first order in the time step") {
  const auto run = [](const mps::Real dt) {
    MixingColumn column({.levels = 20,
                         .surface_temperature_k = 315.0,
                         .surface_heat_capacity_j_m2_k = 2e6,
                         .stellar_flux_w_m2 = 1361.0,
                         .radiation = true});
    const int steps = static_cast<int>(std::lround(3600.0 / dt));
    mps::Real time_s = 0.0;
    for (int step = 0; step < steps; ++step) {
      column.step(time_s, dt);
      time_s += dt;
    }
    std::vector<mps::Real> result = column.potential_temperature_k;
    result.push_back(column.surface_temperature_k);
    return result;
  };
  const auto reference = run(3600.0 / 256.0);
  const auto error = [&](const mps::Real dt) {
    const auto values = run(dt);
    mps::Real result = 0.0;
    for (std::size_t index = 0; index < values.size(); ++index)
      result = std::max(result, std::abs(values[index] - reference[index]));
    return result;
  };
  const mps::Real coarse = error(3600.0 / 8.0);
  const mps::Real medium = error(3600.0 / 16.0);
  const mps::Real fine = error(3600.0 / 32.0);
  MPS_CHECK(coarse > 0.0 && medium > 0.0 && fine > 0.0);
  MPS_CHECK(std::log2(coarse / medium) >= 0.8);
  MPS_CHECK(std::log2(medium / fine) >= 0.8);
}

MPS_TEST_CASE("vertical refinement converges the near-surface mixing response") {
  const auto run = [](const std::size_t levels) {
    MixingColumn column({.levels = levels,
                         .surface_temperature_k = 320.0,
                         .surface_heat_capacity_j_m2_k = 1e9});
    mps::Real time_s = 0.0;
    for (int step = 0; step < 120; ++step) {
      column.step(time_s, 60.0);
      time_s += 60.0;
    }
    return column;
  };
  const auto coarse = run(20);
  const auto medium = run(40);
  const auto fine = run(80);
  const mps::Real coarse_error =
      std::abs(coarse.bulk.diagnostics.boundary_layer_height_m -
               fine.bulk.diagnostics.boundary_layer_height_m);
  const mps::Real medium_error =
      std::abs(medium.bulk.diagnostics.boundary_layer_height_m -
               fine.bulk.diagnostics.boundary_layer_height_m);
  MPS_CHECK(medium_error < coarse_error);
  MPS_CHECK(medium_error < 0.1 * fine.bulk.diagnostics.boundary_layer_height_m);
  // The heated mixed layer itself is resolved consistently across the refinement.
  MPS_CHECK(std::abs(medium.mixed_depth_m() - fine.mixed_depth_m()) <
            std::abs(coarse.mixed_depth_m() - fine.mixed_depth_m()) + 1e-9);
  MPS_CHECK(fine.bulk.diagnostics.boundary_layer_height_m > 0.0);
}

int main() { return mps::test::run_all(); }
