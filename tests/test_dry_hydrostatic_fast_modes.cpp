#include <algorithm>
#include <cmath>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "myplanetsim/dynamics/dry_hydrostatic_fast_modes.hpp"
#include "myplanetsim/dynamics/dry_hydrostatic_fast_operator.hpp"
#include "support/test.hpp"

namespace {

constexpr mps::PlanetParameters planet{6371220.0, 7.29212e-5, 9.80616,
                                       287.0,     1004.0,     100000.0};

[[nodiscard]] mps::AtmosphericHybridCoordinate coordinate() {
  return mps::AtmosphericHybridCoordinate(
      {.a_half_pa = {1000.0, 750.0, 500.0, 250.0, 0.0},
       .b_half = {0.0, 0.25, 0.5, 0.75, 1.0}},
      80000.0, 120000.0, 100.0);
}

constexpr mps::SemiImplicitParameters parameters{
    .reference_surface_pressure_pa = 100000.0,
    .reference_temperature_k = 280.0,
    .implicit_weight = 0.5,
    .wave_cfl_threshold = 0.45,
    .maximum_implicit_modes = 4,
    .nonlinear_relative_tolerance = 1.0e-8,
    .nonlinear_maximum_iterations = 4,
    .linear_relative_tolerance = 1.0e-8,
    .linear_absolute_tolerance = 1.0e-12,
    .linear_maximum_iterations = 40,
    .gmres_restart = 20,
    .minimum_time_step_s = 60.0};

}  // namespace

MPS_TEST_CASE("reference column is isothermal and hydrostatic") {
  const auto reference =
      mps::make_dry_hydrostatic_reference_column(coordinate(), planet, parameters);
  MPS_CHECK_EQ(reference.geometry.air_mass_kg_m2.size(), 4U);
  mps::Real pressure_thickness = 0.0;
  for (std::size_t level = 0; level < 4; ++level) {
    MPS_CHECK_NEAR(
        reference.potential_temperature_k[level] * reference.geometry.exner_full[level],
        parameters.reference_temperature_k, 2.0e-13);
    MPS_CHECK_NEAR(reference.potential_temperature_mass_k_kg_m2[level],
                   reference.geometry.air_mass_kg_m2[level] *
                       reference.potential_temperature_k[level],
                   1.0e-12);
    MPS_CHECK_NEAR(reference.hydrostatic.residual_m2_s2[level], 0.0, 1.0e-10);
    pressure_thickness += reference.geometry.delta_pressure_pa[level];
  }
  MPS_CHECK_NEAR(pressure_thickness,
                 parameters.reference_surface_pressure_pa -
                     reference.geometry.pressure_half_pa.front(),
                 1.0e-10);
}

MPS_TEST_CASE("external mode has the registered Lamb speed and round trips") {
  const auto reference =
      mps::make_dry_hydrostatic_reference_column(coordinate(), planet, parameters);
  const auto modes = mps::make_dry_hydrostatic_external_mode(reference, planet);
  MPS_CHECK_EQ(modes.mode_count(), 1U);
  const mps::Real gamma = planet.heat_capacity_cp_j_kg_k /
                          (planet.heat_capacity_cp_j_kg_k - planet.gas_constant_j_kg_k);
  MPS_CHECK_NEAR(modes.phase_speed_m_s[0],
                 std::sqrt(gamma * planet.gas_constant_j_kg_k *
                           parameters.reference_temperature_k),
                 1.0e-13);
  const std::vector<mps::Real> external_profile(4, 3.0);
  std::vector<mps::Real> modal(1);
  std::vector<mps::Real> reconstructed(4);
  mps::project_onto_vertical_modes(modes, external_profile, modal);
  mps::reconstruct_from_vertical_modes(modes, modal, reconstructed);
  for (const auto value : reconstructed) MPS_CHECK_NEAR(value, 3.0, 1.0e-14);
}

MPS_TEST_CASE("mode selection uses the cell Courant sum and enforces its cap") {
  const auto reference =
      mps::make_dry_hydrostatic_reference_column(coordinate(), planet, parameters);
  const auto modes = mps::make_dry_hydrostatic_external_mode(reference, planet);
  const mps::CubedSphereGrid grid(12, planet.radius_m);
  MPS_CHECK(mps::select_implicit_vertical_modes(grid, modes, 100.0, 0.45, 1).empty());
  const auto selected =
      mps::select_implicit_vertical_modes(grid, modes, 1800.0, 0.45, 1);
  MPS_CHECK_EQ(selected.size(), 1U);
  MPS_CHECK_EQ(selected.front(), 0U);
  MPS_CHECK_THROWS_AS(mps::select_implicit_vertical_modes(grid, modes, 1800.0, 0.45, 0),
                      std::invalid_argument);
}

MPS_TEST_CASE("full vertical modes diagonalize the reference fast structure") {
  const auto c = coordinate();
  const auto reference =
      mps::make_dry_hydrostatic_reference_column(c, planet, parameters);
  const auto fast = mps::make_dry_hydrostatic_fast_operator(c, planet, reference);
  const auto modes = mps::make_dry_hydrostatic_vertical_modes(reference, planet, fast);
  MPS_CHECK_EQ(modes.mode_count(), 4U);
  MPS_CHECK_EQ(modes.vertical_structure_m2_s2.size(), 16U);
  for (std::size_t mode = 0; mode < modes.mode_count(); ++mode) {
    MPS_CHECK(modes.phase_speed_m_s[mode] > 0.0);
    if (mode > 0)
      MPS_CHECK(modes.phase_speed_m_s[mode - 1] >= modes.phase_speed_m_s[mode]);
    const auto eigenvector = modes.eigenvector(mode);
    mps::Real residual = 0.0;
    mps::Real scale = 0.0;
    for (std::size_t row = 0; row < modes.levels; ++row) {
      mps::Real actual = 0.0;
      for (std::size_t column = 0; column < modes.levels; ++column)
        actual += modes.vertical_structure_m2_s2[row * modes.levels + column] *
                  eigenvector[column];
      const mps::Real expected = modes.eigenvalue_m2_s2[mode] * eigenvector[row];
      residual = std::max(residual, std::abs(actual - expected));
      scale = std::max(scale, std::abs(expected));
    }
    MPS_CHECK(residual < 2.0e-8 * std::max(1.0, scale));
  }
  const std::vector<mps::Real> profile = {1.0, -2.0, 0.5, 4.0};
  std::vector<mps::Real> modal(4);
  std::vector<mps::Real> reconstructed(4);
  mps::project_onto_vertical_modes(modes, profile, modal);
  mps::reconstruct_from_vertical_modes(modes, modal, reconstructed);
  for (std::size_t level = 0; level < profile.size(); ++level)
    MPS_CHECK_NEAR(reconstructed[level], profile[level], 2.0e-11);
  std::ostringstream metadata;
  mps::write_dry_hydrostatic_vertical_mode_metadata(metadata, modes);
  MPS_CHECK(metadata.str().find("semi_implicit.vertical_mode_count = 4\n") !=
            std::string::npos);
  MPS_CHECK(metadata.str().find("semi_implicit.mode_0_phase_speed_m_s = ") !=
            std::string::npos);
}

int main() { return mps::test::run_all(); }
