#include "myplanetsim/dynamics/surface_boundary.hpp"

#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <locale>
#include <map>
#include <numbers>
#include <ranges>
#include <sstream>
#include <stdexcept>
#include <utility>

#include "myplanetsim/core/validation.hpp"

namespace mps {
namespace {

struct EarthSourceGrid {
  std::vector<Real> longitudes_deg;
  std::vector<Real> latitudes_deg;
  std::vector<Real> mean_height_m;
  std::vector<Real> land_fraction;
};

[[nodiscard]] Real parse_number(const std::string_view text) {
  Real value = 0.0;
  const auto parsed = std::from_chars(text.data(), text.data() + text.size(), value);
  if (parsed.ec != std::errc{} || parsed.ptr != text.data() + text.size() ||
      !std::isfinite(value))
    throw std::runtime_error("Earth surface CSV contains an invalid number");
  return value;
}

[[nodiscard]] EarthSourceGrid parse_earth_source(const std::string& bytes) {
  std::istringstream input(bytes);
  input.imbue(std::locale::classic());
  std::string line;
  if (!std::getline(input, line) ||
      line != "longitude_deg,latitude_deg,mean_surface_height_m,land_fraction")
    throw std::runtime_error("Earth surface CSV header is invalid");
  std::vector<Real> longitude_order;
  std::vector<Real> latitude_order;
  std::map<std::pair<Real, Real>, std::pair<Real, Real>> records;
  while (std::getline(input, line)) {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (line.empty()) throw std::runtime_error("Earth surface CSV has an empty row");
    std::array<std::string_view, 4> fields{};
    std::size_t begin = 0;
    for (std::size_t field = 0; field < fields.size(); ++field) {
      const auto comma = line.find(',', begin);
      if ((field + 1 < fields.size()) != (comma != std::string::npos))
        throw std::runtime_error("Earth surface CSV row must have four fields");
      fields[field] = std::string_view(line).substr(
          begin, comma == std::string::npos ? comma : comma - begin);
      begin = comma == std::string::npos ? line.size() : comma + 1;
    }
    const Real longitude = parse_number(fields[0]);
    const Real latitude = parse_number(fields[1]);
    const Real height = parse_number(fields[2]);
    const Real fraction = parse_number(fields[3]);
    if (latitude < -90.0 || latitude > 90.0 || fraction < 0.0 || fraction > 1.0 ||
        (fraction == 0.0 && height != 0.0))
      throw std::runtime_error("Earth surface CSV value is out of range");
    if (!records.emplace(std::pair{longitude, latitude}, std::pair{height, fraction})
             .second)
      throw std::runtime_error("Earth surface CSV contains a duplicate cell");
    if (std::ranges::find(longitude_order, longitude) == longitude_order.end())
      longitude_order.push_back(longitude);
    if (std::ranges::find(latitude_order, latitude) == latitude_order.end())
      latitude_order.push_back(latitude);
  }
  if (input.bad() || longitude_order.size() < 2 || latitude_order.size() < 2 ||
      !std::ranges::is_sorted(longitude_order) ||
      !std::ranges::is_sorted(latitude_order) ||
      records.size() != longitude_order.size() * latitude_order.size())
    throw std::runtime_error("Earth surface CSV is not a complete sorted grid");
  EarthSourceGrid result{longitude_order, latitude_order,
                         std::vector<Real>(records.size()),
                         std::vector<Real>(records.size())};
  for (std::size_t longitude = 0; longitude < longitude_order.size(); ++longitude) {
    for (std::size_t latitude = 0; latitude < latitude_order.size(); ++latitude) {
      const auto [height, fraction] =
          records.at({longitude_order[longitude], latitude_order[latitude]});
      const auto index = longitude * latitude_order.size() + latitude;
      result.mean_height_m[index] = height;
      result.land_fraction[index] = fraction;
    }
  }
  return result;
}

[[nodiscard]] std::pair<Real, Real> sample_earth(const EarthSourceGrid& source,
                                                 const Vec3 position) {
  Real longitude =
      std::atan2(position.y, position.x) * 180.0 / std::numbers::pi_v<Real>;
  while (longitude < source.longitudes_deg.front()) longitude += 360.0;
  while (longitude >= source.longitudes_deg.front() + 360.0) longitude -= 360.0;
  const Real latitude =
      std::asin(std::clamp(position.z, -1.0, 1.0)) * 180.0 / std::numbers::pi_v<Real>;
  const auto upper_longitude =
      std::ranges::upper_bound(source.longitudes_deg, longitude);
  const std::size_t longitude1 =
      upper_longitude == source.longitudes_deg.end()
          ? 0
          : static_cast<std::size_t>(upper_longitude - source.longitudes_deg.begin());
  const std::size_t longitude0 =
      longitude1 == 0 ? source.longitudes_deg.size() - 1 : longitude1 - 1;
  const Real left = source.longitudes_deg[longitude0] - (longitude1 == 0 ? 360.0 : 0.0);
  const Real right = source.longitudes_deg[longitude1];
  const Real adjusted_longitude =
      longitude1 == 0 && longitude > right ? longitude - 360.0 : longitude;
  const Real longitude_weight = (adjusted_longitude - left) / (right - left);

  const auto upper_latitude = std::ranges::lower_bound(source.latitudes_deg, latitude);
  std::size_t latitude1 =
      static_cast<std::size_t>(upper_latitude - source.latitudes_deg.begin());
  if (latitude1 == source.latitudes_deg.size()) latitude1 -= 1;
  const std::size_t latitude0 = latitude1 == 0 ? 0 : latitude1 - 1;
  const Real latitude_weight =
      latitude0 == latitude1
          ? 0.0
          : std::clamp(
                (latitude - source.latitudes_deg[latitude0]) /
                    (source.latitudes_deg[latitude1] - source.latitudes_deg[latitude0]),
                0.0, 1.0);
  const auto interpolate = [&](const std::vector<Real>& field) {
    const auto at = [&](const std::size_t i, const std::size_t j) {
      return field[i * source.latitudes_deg.size() + j];
    };
    const Real lower = (1.0 - longitude_weight) * at(longitude0, latitude0) +
                       longitude_weight * at(longitude1, latitude0);
    const Real upper = (1.0 - longitude_weight) * at(longitude0, latitude1) +
                       longitude_weight * at(longitude1, latitude1);
    return (1.0 - latitude_weight) * lower + latitude_weight * upper;
  };
  return {interpolate(source.mean_height_m),
          std::clamp(interpolate(source.land_fraction), 0.0, 1.0)};
}

}  // namespace

SurfaceBoundary::SurfaceBoundary(std::vector<Real> surface_geopotential_m2_s2,
                                 std::vector<Real> land_fraction, std::string source_id,
                                 std::string source_fingerprint)
    : surface_geopotential_m2_s2_(std::move(surface_geopotential_m2_s2)),
      land_fraction_(std::move(land_fraction)),
      source_id_(std::move(source_id)),
      source_fingerprint_(std::move(source_fingerprint)) {
  if (surface_geopotential_m2_s2_.empty() ||
      surface_geopotential_m2_s2_.size() != land_fraction_.size()) {
    throw std::invalid_argument("surface boundary fields must have equal nonzero size");
  }
  for (const Real value : surface_geopotential_m2_s2_) {
    require_finite(value, "surface geopotential");
  }
  for (const Real fraction : land_fraction_) {
    require_finite(fraction, "land fraction");
    if (fraction < 0.0 || fraction > 1.0) {
      throw std::invalid_argument("land fraction must be in [0, 1]");
    }
  }
  if (source_id_.empty() || source_fingerprint_.empty()) {
    throw std::invalid_argument("surface boundary provenance must not be empty");
  }
}

SurfaceBoundary make_surface_boundary(const SurfaceParameters& surface,
                                      const OrographyParameters& orography,
                                      const CubedSphereGrid& grid,
                                      const Real gravity_m_s2,
                                      const std::filesystem::path& source_directory) {
  require_finite(surface.uniform_land_fraction, "surface.uniform_land_fraction");
  if (surface.geography == SurfaceGeography::kEarth) {
    if (orography.kind != OrographyKind::kFlat)
      throw std::invalid_argument("Earth geography owns orography");
    if (surface.quadrature_order <= 0 || surface.quadrature_order > 8 ||
        surface.smoothing_passes < 0)
      throw std::invalid_argument("Earth surface remap options are invalid");
    const auto path = source_directory / surface.input_file;
    std::ifstream input(path, std::ios::binary);
    if (!input)
      throw std::runtime_error("unable to open Earth surface CSV: " + path.string());
    const std::string bytes((std::istreambuf_iterator<char>(input)),
                            std::istreambuf_iterator<char>());
    if (input.bad()) throw std::runtime_error("failed while reading Earth surface CSV");
    const std::string fingerprint = fnv1a64_hex(bytes);
    if (fingerprint != surface.input_fingerprint_fnv1a64)
      throw std::runtime_error("Earth surface CSV fingerprint mismatch");
    const auto source = parse_earth_source(bytes);
    std::vector<Real> geopotential(grid.cell_count());
    std::vector<Real> land_fraction(grid.cell_count());
    for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
      const auto& geometry = grid.cells()[cell];
      std::array<Vec3, 4> corners{};
      for (std::size_t corner = 0; corner < corners.size(); ++corner)
        corners[corner] = grid.vertices()[geometry.vertices[corner]].position;
      Real height_sum = 0.0;
      Real fraction_sum = 0.0;
      for (Index j = 0; j < surface.quadrature_order; ++j) {
        for (Index i = 0; i < surface.quadrature_order; ++i) {
          const Real u = (static_cast<Real>(i) + 0.5) /
                         static_cast<Real>(surface.quadrature_order);
          const Real v = (static_cast<Real>(j) + 0.5) /
                         static_cast<Real>(surface.quadrature_order);
          const Vec3 point = normalize((1.0 - u) * (1.0 - v) * corners[0] +
                                       u * (1.0 - v) * corners[1] + u * v * corners[2] +
                                       (1.0 - u) * v * corners[3]);
          const auto [height, fraction] = sample_earth(source, point);
          height_sum += height;
          fraction_sum += fraction;
        }
      }
      const Real count =
          static_cast<Real>(surface.quadrature_order * surface.quadrature_order);
      land_fraction[cell] = std::clamp(fraction_sum / count, 0.0, 1.0);
      geopotential[cell] = gravity_m_s2 * height_sum / count;
      if (land_fraction[cell] == 0.0) geopotential[cell] = 0.0;
    }
    geopotential =
        smooth_surface_geopotential(grid, geopotential, surface.smoothing_passes);
    return SurfaceBoundary(std::move(geopotential), std::move(land_fraction),
                           "earth:etopo2022-gshhg2.3.7", fingerprint);
  }
  if (surface.uniform_land_fraction < 0.0 || surface.uniform_land_fraction > 1.0) {
    throw std::invalid_argument("surface.uniform_land_fraction must be in [0, 1]");
  }
  const auto terrain =
      make_surface_orography(orography, grid, gravity_m_s2, source_directory);
  std::ostringstream identity;
  identity.imbue(std::locale::classic());
  identity << std::setprecision(std::numeric_limits<Real>::max_digits10)
           << "uniform:" << surface.uniform_land_fraction << ':'
           << terrain.source_fingerprint();
  const std::string source_id = "uniform+" + terrain.source_fingerprint();
  return SurfaceBoundary(
      std::vector<Real>(terrain.surface_geopotential_m2_s2().begin(),
                        terrain.surface_geopotential_m2_s2().end()),
      std::vector<Real>(grid.cell_count(), surface.uniform_land_fraction), source_id,
      fnv1a64_hex(identity.str()));
}

Real mixed_surface_heat_capacity(const Real land_fraction,
                                 const Real land_heat_capacity_j_m2_k,
                                 const Real ocean_heat_capacity_j_m2_k) {
  require_finite(land_fraction, "land fraction");
  if (land_fraction < 0.0 || land_fraction > 1.0) {
    throw std::invalid_argument("land fraction must be in [0, 1]");
  }
  require_positive(land_heat_capacity_j_m2_k, "land surface heat capacity");
  require_positive(ocean_heat_capacity_j_m2_k, "ocean surface heat capacity");
  return land_fraction * land_heat_capacity_j_m2_k +
         (1.0 - land_fraction) * ocean_heat_capacity_j_m2_k;
}

}  // namespace mps
