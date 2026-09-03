#pragma once

#include "myplanetsim/config/experiment_config.hpp"
#include "myplanetsim/dynamics/shallow_water_state.hpp"

namespace mps {

[[nodiscard]] Real linear_wave_frequency(Real radius_m, Real gravity_m_s2,
                                         Real mean_depth_m);
[[nodiscard]] ShallowWaterState make_linear_wave_state(
    const CubedSphereGrid& grid, Real time_s, Real gravity_m_s2, Real mean_depth_m,
    Real relative_amplitude = 1.0e-4);
[[nodiscard]] ShallowWaterState make_geostrophic_adjustment_state(
    const CubedSphereGrid& grid, Real time_s, Real mean_depth_m,
    Real relative_amplitude = 1.0e-2);
[[nodiscard]] ShallowWaterState make_williamson2_state(const CubedSphereGrid& grid,
                                                       Real time_s, Real gravity_m_s2,
                                                       Real rotation_rate_rad_s,
                                                       Real reference_depth_m,
                                                       Real maximum_velocity_m_s,
                                                       Vec3 flow_axis);
[[nodiscard]] ShallowWaterState make_williamson6_state(const CubedSphereGrid& grid,
                                                       Real time_s, Real gravity_m_s2,
                                                       Real rotation_rate_rad_s,
                                                       Real reference_depth_m,
                                                       Real velocity_scale_m_s,
                                                       Vec3 flow_axis);
[[nodiscard]] ShallowWaterState make_galewsky_state(
    const CubedSphereGrid& grid, Real time_s, Real gravity_m_s2,
    Real rotation_rate_rad_s, Real mean_depth_m, Real maximum_velocity_m_s,
    Vec3 flow_axis, bool add_height_perturbation = true);

struct SphericalWaveMode {
  Real amplitude;
  Real phase_rad;
};

[[nodiscard]] SphericalWaveMode diagnose_depth_wave_mode(const CubedSphereGrid& grid,
                                                         const ShallowWaterState& state,
                                                         Vec3 axis, int wavenumber);
[[nodiscard]] ShallowWaterState make_shallow_water_initial_state(
    const CubedSphereGrid& grid, const ExperimentConfig& config);

}  // namespace mps
