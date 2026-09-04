#include "myplanetsim/dynamics/surface_orography.hpp"

#include <cmath>
#include <algorithm>
#include <charconv>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <locale>
#include <map>
#include <numbers>
#include <set>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace mps {
namespace {

struct LatLonTerrain {
  std::vector<Real> longitudes_deg;
  std::vector<Real> latitudes_deg;
  std::vector<Real> height_m;
};

Real parse_csv_real(const std::string_view text) {
  Real value = 0;
  const auto parsed = std::from_chars(text.data(), text.data() + text.size(), value);
  if (parsed.ec != std::errc{} || parsed.ptr != text.data() + text.size() ||
      !std::isfinite(value))
    throw std::runtime_error("terrain CSV contains an invalid number");
  return value;
}

LatLonTerrain parse_latlon_csv(const std::string& bytes) {
  std::istringstream input(bytes);
  input.imbue(std::locale::classic());
  std::string line;
  if (!std::getline(input, line))
    throw std::runtime_error("terrain CSV is empty");
  if (!line.empty() && line.back() == '\r') line.pop_back();
  if (line != "longitude_deg,latitude_deg,height_m")
    throw std::runtime_error("terrain CSV header is invalid");
  std::vector<Real> longitude_order;
  std::vector<Real> latitude_order;
  std::map<std::pair<Real, Real>, Real> records;
  while (std::getline(input, line)) {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (line.empty()) throw std::runtime_error("terrain CSV contains an empty row");
    const auto first = line.find(',');
    const auto second = first == std::string::npos ? first : line.find(',', first + 1);
    if (first == std::string::npos || second == std::string::npos ||
        line.find(',', second + 1) != std::string::npos)
      throw std::runtime_error("terrain CSV row must contain exactly three fields");
    const Real longitude = parse_csv_real(
        std::string_view(line).substr(0, first));
    const Real latitude = parse_csv_real(
        std::string_view(line).substr(first + 1, second - first - 1));
    const Real height = parse_csv_real(std::string_view(line).substr(second + 1));
    if (latitude < -90 || latitude > 90 || height < 0)
      throw std::runtime_error("terrain CSV latitude or height is out of range");
    if (!records.emplace(std::pair{longitude, latitude}, height).second)
      throw std::runtime_error("terrain CSV contains a duplicate grid point");
    if (std::ranges::find(longitude_order, longitude) == longitude_order.end())
      longitude_order.push_back(longitude);
    if (std::ranges::find(latitude_order, latitude) == latitude_order.end())
      latitude_order.push_back(latitude);
  }
  if (input.bad() || longitude_order.size() < 2 || latitude_order.size() < 2)
    throw std::runtime_error("terrain CSV grid is incomplete");
  if (!std::ranges::is_sorted(longitude_order, std::less{}) ||
      std::ranges::adjacent_find(longitude_order) != longitude_order.end() ||
      !std::ranges::is_sorted(latitude_order, std::less{}) ||
      std::ranges::adjacent_find(latitude_order) != latitude_order.end())
    throw std::runtime_error("terrain CSV axes must be strictly increasing");
  if (!(longitude_order.back() < longitude_order.front() + 360.0))
    throw std::runtime_error("terrain CSV longitude must contain one periodic cycle");
  if (records.size() != longitude_order.size() * latitude_order.size())
    throw std::runtime_error("terrain CSV tensor product has missing points");
  LatLonTerrain result{longitude_order, latitude_order,
                       std::vector<Real>(records.size())};
  for (std::size_t longitude = 0; longitude < longitude_order.size(); ++longitude)
    for (std::size_t latitude = 0; latitude < latitude_order.size(); ++latitude)
      result.height_m[longitude * latitude_order.size() + latitude] =
          records.at({longitude_order[longitude], latitude_order[latitude]});
  return result;
}

std::vector<Real> interpolate_latlon(const LatLonTerrain& terrain,
                                     const CubedSphereGrid& grid,
                                     const Real gravity_m_s2) {
  std::vector<Real> result(grid.cell_count());
  const auto& longitude = terrain.longitudes_deg;
  const auto& latitude = terrain.latitudes_deg;
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
    const Vec3 position = grid.cells()[cell].center;
    Real lon = std::atan2(position.y, position.x) * 180.0 /
               std::numbers::pi_v<Real>;
    while (lon < longitude.front()) lon += 360.0;
    while (lon >= longitude.front() + 360.0) lon -= 360.0;
    const Real lat = std::asin(std::clamp(position.z, -1.0, 1.0)) * 180.0 /
                     std::numbers::pi_v<Real>;
    if (lat < latitude.front() || lat > latitude.back())
      throw std::runtime_error("terrain CSV latitude does not cover the model grid");
    auto upper_lon = std::ranges::upper_bound(longitude, lon);
    std::size_t lon0 = 0;
    std::size_t lon1 = 0;
    Real lon_left = 0;
    Real lon_right = 0;
    if (upper_lon == longitude.end()) {
      lon0 = longitude.size() - 1;
      lon1 = 0;
      lon_left = longitude[lon0];
      lon_right = longitude.front() + 360.0;
    } else {
      lon1 = static_cast<std::size_t>(upper_lon - longitude.begin());
      lon0 = lon1 == 0 ? longitude.size() - 1 : lon1 - 1;
      lon_left = lon1 == 0 ? longitude[lon0] - 360.0 : longitude[lon0];
      lon_right = longitude[lon1];
    }
    auto upper_lat = std::ranges::lower_bound(latitude, lat);
    std::size_t lat1 = static_cast<std::size_t>(upper_lat - latitude.begin());
    if (lat1 == latitude.size()) lat1 = latitude.size() - 1;
    const std::size_t lat0 = lat1 == 0 ? 0 : lat1 - 1;
    const Real wx = (lon - lon_left) / (lon_right - lon_left);
    const Real wy = lat0 == lat1 ? 0.0 :
        (lat - latitude[lat0]) / (latitude[lat1] - latitude[lat0]);
    const auto sample = [&](const std::size_t i, const std::size_t j) {
      return terrain.height_m[i * latitude.size() + j];
    };
    const Real lower = (1 - wx) * sample(lon0, lat0) + wx * sample(lon1, lat0);
    const Real upper = (1 - wx) * sample(lon0, lat1) + wx * sample(lon1, lat1);
    result[cell] = gravity_m_s2 * ((1 - wy) * lower + wy * upper);
  }
  return result;
}

}  // namespace

std::string fnv1a64_hex(const std::string_view bytes) {
  std::uint64_t hash = 14695981039346656037ULL;
  for (const unsigned char byte : bytes) {
    hash ^= byte;
    hash *= 1099511628211ULL;
  }
  std::ostringstream output;
  output << std::hex << std::setfill('0') << std::setw(16) << hash;
  return output.str();
}

std::vector<Real> smooth_surface_geopotential(
    const CubedSphereGrid& grid, const std::span<const Real> values,
    const Index passes) {
  if (values.size() != grid.cell_count() || passes < 0)
    throw std::invalid_argument("terrain smoothing arguments are invalid");
  std::vector<Real> current(values.begin(), values.end());
  for (const auto value : current)
    if (!std::isfinite(value))
      throw std::invalid_argument("terrain smoothing field must be finite");
  for (Index pass = 0; pass < passes; ++pass) {
    std::vector<Real> integrated_increment(grid.cell_count());
    for (const auto& edge : grid.edges()) {
      const auto left = grid.cell_index(edge.left_cell);
      const auto right = grid.cell_index(edge.right_cell);
      const Real capacity = 0.5 * std::min(
          grid.cells()[left].area_m2 /
              static_cast<Real>(grid.cell_edges(edge.left_cell).size()),
          grid.cells()[right].area_m2 /
              static_cast<Real>(grid.cell_edges(edge.right_cell).size()));
      const Real exchange = capacity * (current[right] - current[left]);
      integrated_increment[left] += exchange;
      integrated_increment[right] -= exchange;
    }
    for (std::size_t cell = 0; cell < current.size(); ++cell)
      current[cell] += integrated_increment[cell] / grid.cells()[cell].area_m2;
  }
  return current;
}

SurfaceOrography::SurfaceOrography(
    std::vector<Real> surface_geopotential_m2_s2,
    std::string source_fingerprint)
    : surface_geopotential_m2_s2_(std::move(surface_geopotential_m2_s2)),
      source_fingerprint_(std::move(source_fingerprint)) {
  if (surface_geopotential_m2_s2_.empty())
    throw std::invalid_argument("surface orography must contain cells");
  for (const auto value : surface_geopotential_m2_s2_)
    if (!std::isfinite(value))
      throw std::invalid_argument("surface orography must be finite");
  if (source_fingerprint_.empty())
    throw std::invalid_argument("surface orography fingerprint must not be empty");
}

Real dcmip_2_0_0_surface_height_m(const Vec3 position) {
  if (!is_finite(position) || !(norm_squared(position) > 0.0))
    throw std::invalid_argument("DCMIP terrain position must be finite and nonzero");
  const Vec3 unit = normalize(position);
  constexpr Vec3 mountain_centre{0.0, -1.0, 0.0};
  constexpr Real height_m = 2000.0;
  constexpr Real radius = 3.0 * std::numbers::pi_v<Real> / 4.0;
  constexpr Real half_width = std::numbers::pi_v<Real> / 16.0;
  const Real distance = safe_angle(unit, mountain_centre);
  if (!(distance < radius)) return 0.0;
  const Real envelope = 0.5 * (1.0 + std::cos(std::numbers::pi_v<Real> *
                                               distance / radius));
  const Real oscillation =
      std::cos(std::numbers::pi_v<Real> * distance / half_width);
  return height_m * envelope * oscillation * oscillation;
}

Real williamson5_surface_height_m(const Vec3 position) {
  if (!is_finite(position) || !(norm_squared(position) > 0.0))
    throw std::invalid_argument(
        "Williamson 5 terrain position must be finite and nonzero");
  const Vec3 unit = normalize(position);
  const Real longitude = std::atan2(unit.y, unit.x);
  const Real latitude = std::asin(std::clamp(unit.z, -1.0, 1.0));
  constexpr Real centre_longitude = -0.5 * std::numbers::pi_v<Real>;
  constexpr Real centre_latitude = std::numbers::pi_v<Real> / 6.0;
  constexpr Real radius = std::numbers::pi_v<Real> / 9.0;
  const Real longitude_distance =
      std::remainder(longitude - centre_longitude, 2.0 * std::numbers::pi_v<Real>);
  const Real distance =
      std::hypot(longitude_distance, latitude - centre_latitude);
  return distance < radius ? 2000.0 * (1.0 - distance / radius) : 0.0;
}

Real linear_bell_surface_height_m(const Vec3 position, const Real peak_height_m) {
  if (!is_finite(position) || !(norm_squared(position) > 0.0) ||
      !(peak_height_m > 0.0) || !std::isfinite(peak_height_m))
    throw std::invalid_argument("linear bell arguments are invalid");
  constexpr Vec3 centre{1.0, 0.0, 0.0};
  constexpr Real radius = std::numbers::pi_v<Real> / 6.0;
  const Real distance = safe_angle(normalize(position), centre);
  if (!(distance < radius)) return 0.0;
  const Real taper = std::cos(0.5 * std::numbers::pi_v<Real> * distance / radius);
  return peak_height_m * taper * taper;
}

Real jw06_surface_geopotential_m2_s2(const Vec3 position) {
  if (!is_finite(position) || !(norm_squared(position) > 0.0))
    throw std::invalid_argument("JW06 terrain position must be finite and nonzero");
  constexpr Real u0 = 35.0;
  constexpr Real eta0 = 0.252;
  constexpr Real radius_m = 6.371229e6;
  constexpr Real rotation_rate_rad_s = 7.29212e-5;
  const Real latitude = std::asin(std::clamp(normalize(position).z, -1.0, 1.0));
  const Real sine = std::sin(latitude);
  const Real cosine = std::cos(latitude);
  const Real eta_v = (1.0 - eta0) * 0.5 * std::numbers::pi_v<Real>;
  const Real vertical = std::pow(std::cos(eta_v), 1.5);
  const Real wind_shape =
      -2.0 * std::pow(sine, 6) * (cosine * cosine + 1.0 / 3.0) + 10.0 / 63.0;
  const Real rotation_shape =
      8.0 / 5.0 * std::pow(cosine, 3) * (sine * sine + 2.0 / 3.0) -
      std::numbers::pi_v<Real> / 4.0;
  return u0 * vertical *
         (u0 * vertical * wind_shape +
          radius_m * rotation_rate_rad_s * rotation_shape);
}

SurfaceOrography make_surface_orography(const OrographyParameters& parameters,
                                        const CubedSphereGrid& grid,
                                        const Real gravity_m_s2,
                                        const std::filesystem::path& source_directory) {
  if (!(gravity_m_s2 > 0.0) || !std::isfinite(gravity_m_s2))
    throw std::invalid_argument("surface orography gravity must be positive");
  if (parameters.kind == OrographyKind::kFlat)
    return SurfaceOrography(std::vector<Real>(grid.cell_count(), 0.0),
                            "analytic:flat");
  if (parameters.kind == OrographyKind::kDcmip200) {
    std::vector<Real> geopotential(grid.cell_count());
    for (std::size_t cell = 0; cell < grid.cell_count(); ++cell)
      geopotential[cell] = gravity_m_s2 *
                           dcmip_2_0_0_surface_height_m(grid.cells()[cell].center);
    return SurfaceOrography(std::move(geopotential), "analytic:dcmip_2_0_0");
  }
  if (parameters.kind == OrographyKind::kWilliamson5) {
    std::vector<Real> geopotential(grid.cell_count());
    for (std::size_t cell = 0; cell < grid.cell_count(); ++cell)
      geopotential[cell] = gravity_m_s2 *
                           williamson5_surface_height_m(grid.cells()[cell].center);
    return SurfaceOrography(std::move(geopotential), "analytic:williamson5");
  }
  if (parameters.kind == OrographyKind::kLinearBell) {
    std::vector<Real> geopotential(grid.cell_count());
    for (std::size_t cell = 0; cell < grid.cell_count(); ++cell)
      geopotential[cell] = gravity_m_s2 *
                           linear_bell_surface_height_m(grid.cells()[cell].center);
    return SurfaceOrography(std::move(geopotential), "analytic:linear_bell");
  }
  if (parameters.kind == OrographyKind::kJw06) {
    std::vector<Real> geopotential(grid.cell_count());
    for (std::size_t cell = 0; cell < grid.cell_count(); ++cell)
      geopotential[cell] =
          jw06_surface_geopotential_m2_s2(grid.cells()[cell].center);
    return SurfaceOrography(std::move(geopotential), "analytic:jw06");
  }
  if (parameters.kind == OrographyKind::kLatLonCsv) {
    const auto path = source_directory / parameters.input_file;
    std::ifstream input(path, std::ios::binary);
    if (!input) throw std::runtime_error("unable to open terrain CSV: " + path.string());
    const std::string bytes((std::istreambuf_iterator<char>(input)),
                            std::istreambuf_iterator<char>());
    if (input.bad()) throw std::runtime_error("failed while reading terrain CSV");
    const auto fingerprint = fnv1a64_hex(bytes);
    if (fingerprint != parameters.input_fingerprint_fnv1a64)
      throw std::runtime_error("terrain CSV fingerprint mismatch");
    auto geopotential =
        interpolate_latlon(parse_latlon_csv(bytes), grid, gravity_m_s2);
    geopotential = smooth_surface_geopotential(
        grid, geopotential, parameters.smoothing_passes);
    return SurfaceOrography(std::move(geopotential), fingerprint);
  }
  throw std::invalid_argument("selected orography initializer is not implemented");
}

}  // namespace mps
