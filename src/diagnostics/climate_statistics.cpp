#include "myplanetsim/diagnostics/climate_statistics.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <ostream>
#include <stdexcept>

namespace mps {
namespace {

constexpr Real kRadiansToDegrees = 180.0 / 3.14159265358979323846;

struct InstantaneousMoments {
  Real area = 0.0;
  Real pressure = 0.0;
  Real temperature = 0.0;
  Real zonal_wind = 0.0;
  Real meridional_wind = 0.0;
  Real zonal_wind_squared = 0.0;
  Real meridional_wind_squared = 0.0;
  Real temperature_squared = 0.0;
  Real zonal_meridional_wind = 0.0;
  Real meridional_wind_temperature = 0.0;
  Real layer_mass = 0.0;
  Real meridional_mass_transport = 0.0;
};

[[nodiscard]] Vec3 eastward_unit(const Vec3 position) {
  const Vec3 east{-position.y, position.x, 0.0};
  const Real magnitude = norm(east);
  return magnitude > 0.0 ? east / magnitude : Vec3{0.0, 1.0, 0.0};
}

[[nodiscard]] std::size_t latitude_bin(const Vec3 position, const std::size_t bins) {
  const Real scaled =
      0.5 * (std::clamp(position.z, -1.0, 1.0) + 1.0) * static_cast<Real>(bins);
  return std::min(static_cast<std::size_t>(scaled), bins - 1);
}

}  // namespace

ClimateStatisticsAccumulator::ClimateStatisticsAccumulator(
    const CubedSphereGrid& grid, const std::size_t levels,
    const Real accumulation_start_s)
    : grid_(grid),
      levels_(levels),
      bins_(2 * static_cast<std::size_t>(grid.cells_per_panel())),
      accumulation_start_s_(accumulation_start_s),
      accumulated_(bins_ * levels_) {
  if (levels_ == 0 || !std::isfinite(accumulation_start_s_) ||
      accumulation_start_s_ < 0.0) {
    throw std::invalid_argument("climate statistics dimensions or start time invalid");
  }
}

void ClimateStatisticsAccumulator::observe(const DryHydrostaticState& state,
                                           const DryHydrostaticDerived& derived) {
  const std::size_t volume = grid_.cell_count() * levels_;
  if (derived.cells != grid_.cell_count() || derived.levels != levels_ ||
      state.surface_pressure_pa.size() != derived.cells ||
      derived.pressure_pa.size() != volume || derived.air_mass_kg_m2.size() != volume ||
      derived.temperature_k.size() != volume || derived.velocity_m_s.size() != volume) {
    throw std::invalid_argument("climate statistics state shape mismatch");
  }
  if (has_previous_time_ && state.time_s < previous_time_s_) {
    throw std::invalid_argument("climate statistics time moved backwards");
  }
  if (!has_previous_time_) {
    previous_time_s_ = state.time_s;
    has_previous_time_ = true;
    return;
  }
  const Real sample_start = std::max(previous_time_s_, accumulation_start_s_);
  const Real weight_s = state.time_s - sample_start;
  previous_time_s_ = state.time_s;
  if (!(weight_s > 0.0)) return;

  std::vector<InstantaneousMoments> moments(bins_ * levels_);
  for (std::size_t cell = 0; cell < derived.cells; ++cell) {
    const auto& geometry = grid_.cells()[cell];
    const auto bin = latitude_bin(geometry.center, bins_);
    const Vec3 east = eastward_unit(geometry.center);
    const Vec3 north = cross(geometry.center, east);
    for (std::size_t level = 0; level < levels_; ++level) {
      const auto offset = dry_hydrostatic_offset(cell, level, levels_);
      auto& value = moments[bin * levels_ + level];
      const Real area = geometry.area_m2;
      const Real zonal = dot(derived.velocity_m_s[offset], east);
      const Real meridional = dot(derived.velocity_m_s[offset], north);
      const Real temperature = derived.temperature_k[offset];
      value.area += area;
      value.pressure += area * derived.pressure_pa[offset];
      value.temperature += area * temperature;
      value.zonal_wind += area * zonal;
      value.meridional_wind += area * meridional;
      value.zonal_wind_squared += area * zonal * zonal;
      value.meridional_wind_squared += area * meridional * meridional;
      value.temperature_squared += area * temperature * temperature;
      value.zonal_meridional_wind += area * zonal * meridional;
      value.meridional_wind_temperature += area * meridional * temperature;
      value.layer_mass += area * derived.air_mass_kg_m2[offset];
      value.meridional_mass_transport +=
          area * derived.air_mass_kg_m2[offset] * meridional;
    }
  }

  for (std::size_t bin = 0; bin < bins_; ++bin) {
    const Real lower_mu =
        -1.0 + 2.0 * static_cast<Real>(bin) / static_cast<Real>(bins_);
    const Real upper_mu =
        -1.0 + 2.0 * static_cast<Real>(bin + 1) / static_cast<Real>(bins_);
    const Real meridional_width =
        grid_.radius_m() * (std::asin(upper_mu) - std::asin(lower_mu));
    for (std::size_t level = 0; level < levels_; ++level) {
      const auto offset = bin * levels_ + level;
      const auto& value = moments[offset];
      if (!(value.area > 0.0)) continue;
      const Real inverse_area = 1.0 / value.area;
      const Real mean_u = value.zonal_wind * inverse_area;
      const Real mean_v = value.meridional_wind * inverse_area;
      const Real mean_t = value.temperature * inverse_area;
      auto& out = accumulated_[offset];
      out.pressure_time_integral += weight_s * value.pressure * inverse_area;
      out.temperature_time_integral += weight_s * mean_t;
      out.zonal_wind_time_integral += weight_s * mean_u;
      out.meridional_wind_time_integral += weight_s * mean_v;
      out.eddy_momentum_flux_time_integral +=
          weight_s * (value.zonal_meridional_wind * inverse_area - mean_u * mean_v);
      out.eddy_heat_flux_time_integral +=
          weight_s *
          (value.meridional_wind_temperature * inverse_area - mean_v * mean_t);
      out.eddy_kinetic_energy_time_integral +=
          0.5 * weight_s *
          (value.zonal_wind_squared * inverse_area - mean_u * mean_u +
           value.meridional_wind_squared * inverse_area - mean_v * mean_v);
      out.temperature_variance_time_integral +=
          weight_s * (value.temperature_squared * inverse_area - mean_t * mean_t);
      out.layer_mass_time_integral += weight_s * value.layer_mass;
      out.meridional_mass_transport_time_integral +=
          weight_s * value.meridional_mass_transport / meridional_width;
      out.accumulated_time_s += weight_s;
    }
  }
}

std::vector<ClimateStatisticsRow> ClimateStatisticsAccumulator::rows() const {
  std::vector<ClimateStatisticsRow> result;
  result.reserve(accumulated_.size());
  for (std::size_t bin = 0; bin < bins_; ++bin) {
    const Real lower_mu =
        -1.0 + 2.0 * static_cast<Real>(bin) / static_cast<Real>(bins_);
    const Real upper_mu =
        -1.0 + 2.0 * static_cast<Real>(bin + 1) / static_cast<Real>(bins_);
    Real streamfunction = 0.0;
    for (std::size_t level = 0; level < levels_; ++level) {
      const auto& value = accumulated_[bin * levels_ + level];
      if (value.accumulated_time_s > 0.0) {
        streamfunction +=
            value.meridional_mass_transport_time_integral / value.accumulated_time_s;
      }
      const Real inverse_time =
          value.accumulated_time_s > 0.0 ? 1.0 / value.accumulated_time_s : 0.0;
      result.push_back(
          {.latitude_lower_deg = std::asin(lower_mu) * kRadiansToDegrees,
           .latitude_upper_deg = std::asin(upper_mu) * kRadiansToDegrees,
           .level = level,
           .mean_pressure_pa = value.pressure_time_integral * inverse_time,
           .mean_temperature_k = value.temperature_time_integral * inverse_time,
           .mean_zonal_wind_m_s = value.zonal_wind_time_integral * inverse_time,
           .mean_meridional_wind_m_s =
               value.meridional_wind_time_integral * inverse_time,
           .eddy_momentum_flux_m2_s2 =
               value.eddy_momentum_flux_time_integral * inverse_time,
           .eddy_heat_flux_k_m_s = value.eddy_heat_flux_time_integral * inverse_time,
           .eddy_kinetic_energy_m2_s2 =
               value.eddy_kinetic_energy_time_integral * inverse_time,
           .temperature_variance_k2 =
               value.temperature_variance_time_integral * inverse_time,
           .mass_streamfunction_kg_s = streamfunction,
           .accumulated_time_s = value.accumulated_time_s});
    }
  }
  return result;
}

void write_climate_statistics_csv(std::ostream& output,
                                  const std::vector<ClimateStatisticsRow>& rows) {
  output << "latitude_lower_deg,latitude_upper_deg,level,mean_pressure_pa,"
            "mean_temperature_k,mean_zonal_wind_m_s,mean_meridional_wind_m_s,"
            "eddy_momentum_flux_m2_s2,eddy_heat_flux_k_m_s,"
            "eddy_kinetic_energy_m2_s2,temperature_variance_k2,"
            "mass_streamfunction_kg_s,accumulated_time_s\n";
  output << std::setprecision(std::numeric_limits<Real>::max_digits10);
  for (const auto& row : rows) {
    output << row.latitude_lower_deg << ',' << row.latitude_upper_deg << ','
           << row.level << ',' << row.mean_pressure_pa << ',' << row.mean_temperature_k
           << ',' << row.mean_zonal_wind_m_s << ',' << row.mean_meridional_wind_m_s
           << ',' << row.eddy_momentum_flux_m2_s2 << ',' << row.eddy_heat_flux_k_m_s
           << ',' << row.eddy_kinetic_energy_m2_s2 << ',' << row.temperature_variance_k2
           << ',' << row.mass_streamfunction_kg_s << ',' << row.accumulated_time_s
           << '\n';
  }
  if (!output) throw std::runtime_error("failed while writing climate statistics CSV");
}

}  // namespace mps
