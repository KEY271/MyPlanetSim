#pragma once

#include <iosfwd>
#include <vector>

#include "myplanetsim/dynamics/dry_hydrostatic_state.hpp"
#include "myplanetsim/grid/cubed_sphere_grid.hpp"

namespace mps {

struct ClimateStatisticsRow {
  Real latitude_lower_deg;
  Real latitude_upper_deg;
  std::size_t level;
  Real mean_pressure_pa;
  Real mean_temperature_k;
  Real mean_zonal_wind_m_s;
  Real mean_meridional_wind_m_s;
  Real eddy_momentum_flux_m2_s2;
  Real eddy_heat_flux_k_m_s;
  Real eddy_kinetic_energy_m2_s2;
  Real temperature_variance_k2;
  Real mass_streamfunction_kg_s;
  Real accumulated_time_s;
};

class ClimateStatisticsAccumulator {
 public:
  ClimateStatisticsAccumulator(const CubedSphereGrid& grid, std::size_t levels,
                               Real accumulation_start_s = 200.0 * 86400.0);

  void observe(const DryHydrostaticState& state, const DryHydrostaticDerived& derived);
  [[nodiscard]] std::vector<ClimateStatisticsRow> rows() const;
  [[nodiscard]] std::size_t latitude_bins() const noexcept { return bins_; }

 private:
  struct AccumulatedBinLevel {
    Real pressure_time_integral = 0.0;
    Real temperature_time_integral = 0.0;
    Real zonal_wind_time_integral = 0.0;
    Real meridional_wind_time_integral = 0.0;
    Real eddy_momentum_flux_time_integral = 0.0;
    Real eddy_heat_flux_time_integral = 0.0;
    Real eddy_kinetic_energy_time_integral = 0.0;
    Real temperature_variance_time_integral = 0.0;
    Real layer_mass_time_integral = 0.0;
    Real meridional_mass_transport_time_integral = 0.0;
    Real accumulated_time_s = 0.0;
  };

  const CubedSphereGrid& grid_;
  std::size_t levels_;
  std::size_t bins_;
  Real accumulation_start_s_;
  Real previous_time_s_ = 0.0;
  bool has_previous_time_ = false;
  std::vector<AccumulatedBinLevel> accumulated_;
};

void write_climate_statistics_csv(std::ostream& output,
                                  const std::vector<ClimateStatisticsRow>& rows);

}  // namespace mps
