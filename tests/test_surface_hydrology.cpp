#include <cmath>
#include <utility>
#include <vector>

#include "myplanetsim/physics/boundary_layer.hpp"
#include "myplanetsim/physics/moist_thermodynamics.hpp"
#include "myplanetsim/physics/surface_hydrology.hpp"
#include "support/test.hpp"

namespace {

struct MoistColumn {
  std::vector<double> theta{285.0, 290.0};
  std::vector<mps::Vec3> velocity{{}, {}};
  std::vector<double> vapor{0.005, 0.005};
  std::vector<double> mass{300.0, 700.0};
  std::vector<double> exner_full{1.0, 1.0};
  std::vector<double> exner_half{1.0, 1.0, 1.0};
  std::vector<double> height{1500.0, 500.0};
  std::vector<double> density_half{1.0, 1.0, 1.0};
  std::vector<double> momentum_k{0.0, 0.0, 0.0};
  std::vector<double> heat_k{0.0, 20.0, 0.0};
  std::vector<double> tracer_k{0.0, 20.0, 0.0};
  double surface_temperature = 300.0;

  [[nodiscard]] mps::BoundaryLayerColumnInput input(const double land_fraction,
                                                    const double water,
                                                    const double capacity,
                                                    const double land_conductance,
                                                    const double ocean_conductance,
                                                    const double time_step) const {
    return {.potential_temperature_k = theta,
            .velocity_m_s = velocity,
            .tracer_mixing_ratio = vapor,
            .air_mass_kg_m2 = mass,
            .exner_full = exner_full,
            .exner_half = exner_half,
            .height_full_m = height,
            .density_half_kg_m3 = density_half,
            .eddy_diffusivity_momentum_m2_s = momentum_k,
            .eddy_diffusivity_heat_m2_s = heat_k,
            .eddy_diffusivity_tracer_m2_s = tracer_k,
            .surface_temperature_k = surface_temperature,
            .surface_exner = 1.0,
            .surface_heat_capacity_j_m2_k = 1.0e7,
            .surface_heat_conductance_w_m2_k = 20.0,
            .surface_drag_conductance_kg_m2_s = 0.0,
            .heat_capacity_cp_j_kg_k = 1004.0,
            .time_step_s = time_step,
            .enable_surface_water_exchange = true,
            .surface_pressure_pa = 100000.0,
            .land_fraction = land_fraction,
            .land_water_kg_m2 = water,
            .bucket_capacity_kg_m2 = capacity,
            .bucket_wet_threshold_fraction = 0.75,
            .surface_water_conductance_land_kg_m2_s = land_conductance,
            .surface_water_conductance_ocean_kg_m2_s = ocean_conductance,
            .moist_thermodynamics = {}};
  }
};

[[nodiscard]] std::pair<mps::BoundaryLayerColumnResult, mps::BoundaryLayerColumnResult>
solve_with_both_root_solvers(mps::BoundaryLayerColumnInput input) {
  input.surface_water_root_solver = mps::SurfaceWaterRootSolver::bisection;
  auto bisection = mps::implicit_boundary_layer_column(input);
  input.surface_water_root_solver = mps::SurfaceWaterRootSolver::safeguarded_newton;
  auto newton = mps::implicit_boundary_layer_column(input);
  return {std::move(bisection), std::move(newton)};
}

}  // namespace

MPS_TEST_CASE("ocean evaporation cools the surface and moistens one implicit column") {
  const MoistColumn column;
  const auto result = mps::implicit_boundary_layer_column(
      column.input(0.0, 0.0, 0.0, 0.0, 0.01, 300.0));
  MPS_CHECK(result.surface_water_flux_kg_m2_s > 0.0);
  MPS_CHECK(result.surface_temperature_k < column.surface_temperature);
  MPS_CHECK(result.tracer_mixing_ratio.back() > column.vapor.back());
  MPS_CHECK_NEAR(result.diagnostics.tracer_mass_change_kg_m2,
                 result.diagnostics.evaporation_kg_m2, 1e-11);
  MPS_CHECK_NEAR(result.diagnostics.water_budget_residual_kg_m2, 0.0, 1e-11);
  MPS_CHECK_NEAR(result.diagnostics.moist_enthalpy_budget_residual_j_m2, 0.0, 1e-6);
}

MPS_TEST_CASE(
    "dew warms a full land bucket and overflows through the external ledger") {
  MoistColumn column;
  column.surface_temperature = 280.0;
  column.theta = {290.0, 290.0};
  column.vapor = {0.015, 0.015};
  const auto result = mps::implicit_boundary_layer_column(
      column.input(1.0, 1.0, 1.0, 0.01, 0.0, 300.0));
  MPS_CHECK(result.surface_water_flux_kg_m2_s < 0.0);
  MPS_CHECK(result.surface_temperature_k > column.surface_temperature);
  MPS_CHECK_EQ(result.land_water_kg_m2, 1.0);
  MPS_CHECK(result.diagnostics.runoff_kg_m2_land > 0.0);
  MPS_CHECK_NEAR(result.diagnostics.external_outflow_kg_m2,
                 result.diagnostics.runoff_kg_m2_land, 1e-14);
  MPS_CHECK_NEAR(result.diagnostics.water_budget_residual_kg_m2, 0.0, 1e-11);
  MPS_CHECK_NEAR(result.diagnostics.moist_enthalpy_budget_residual_j_m2, 0.0, 1e-6);
}

MPS_TEST_CASE("empty and nearly empty land buckets enforce the supply constraint") {
  const MoistColumn column;
  const auto empty = mps::implicit_boundary_layer_column(
      column.input(1.0, 0.0, 150.0, 0.02, 0.0, 300.0));
  MPS_CHECK_EQ(empty.surface_water_flux_kg_m2_s, 0.0);
  MPS_CHECK_EQ(empty.land_water_kg_m2, 0.0);

  constexpr double initial_water = 1e-3;
  constexpr double time_step = 300.0;
  const auto limited = mps::implicit_boundary_layer_column(
      column.input(1.0, initial_water, 150.0, 10.0, 0.0, time_step));
  MPS_CHECK(limited.land_water_flux_kg_m2_s >= 0.0);
  MPS_CHECK(time_step * limited.land_water_flux_kg_m2_s <= initial_water);
  MPS_CHECK(limited.land_water_kg_m2 >= 0.0);
  MPS_CHECK_NEAR(limited.diagnostics.water_budget_residual_kg_m2, 0.0, 1e-11);
}

MPS_TEST_CASE("safeguarded surface-water Newton matches bisection") {
  MoistColumn evaporation;
  MoistColumn dew;
  dew.surface_temperature = 280.0;
  dew.theta = {290.0, 290.0};
  dew.vapor = {0.015, 0.015};
  for (auto input : {evaporation.input(0.0, 0.0, 0.0, 0.0, 0.01, 900.0),
                     evaporation.input(1.0, 1e-3, 150.0, 10.0, 0.0, 900.0),
                     dew.input(1.0, 1.0, 1.0, 0.01, 0.0, 900.0)}) {
    const auto [bisection, newton] = solve_with_both_root_solvers(input);
    MPS_CHECK_NEAR(newton.surface_water_flux_kg_m2_s,
                   bisection.surface_water_flux_kg_m2_s, 1e-12);
    MPS_CHECK_NEAR(newton.land_water_kg_m2, bisection.land_water_kg_m2, 1e-9);
    MPS_CHECK_NEAR(newton.surface_temperature_k, bisection.surface_temperature_k, 1e-9);
    MPS_CHECK_NEAR(newton.diagnostics.water_budget_residual_kg_m2, 0.0, 1e-11);
    MPS_CHECK(newton.diagnostics.surface_water_iterations <=
              bisection.diagnostics.surface_water_iterations);
  }
}

MPS_TEST_CASE("mixed land ocean precipitation fills the bucket then routes runoff") {
  const auto update = mps::route_surface_precipitation(0.25, 10.0, 145.0, 150.0);
  MPS_CHECK_EQ(update.final_land_water_kg_m2, 150.0);
  MPS_CHECK_EQ(update.runoff_kg_m2_land, 5.0);
  MPS_CHECK_NEAR(update.ocean_gain_kg_m2_cell, 8.75, 1e-14);
  MPS_CHECK_NEAR(update.water_budget_residual_kg_m2, 0.0, 1e-14);

  const auto ocean = mps::route_surface_precipitation(0.0, 3.0, 0.0, 0.0);
  MPS_CHECK_EQ(ocean.final_land_water_kg_m2, 0.0);
  MPS_CHECK_EQ(ocean.ocean_gain_kg_m2_cell, 3.0);
}

int main() { return mps::test::run_all(); }
