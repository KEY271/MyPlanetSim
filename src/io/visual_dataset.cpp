#include "myplanetsim/io/visual_dataset.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <locale>
#include <sstream>
#include <stdexcept>
#include <system_error>

#include "myplanetsim/core/validation.hpp"
#include "myplanetsim/dynamics/tracer_registry.hpp"

namespace mps {
namespace {

constexpr std::array<char, 8> kTerrainMagic{'M', 'P', 'S', 'T', 'E', 'R', 'R', '1'};
constexpr std::array<char, 8> kMeanMagic{'M', 'P', 'S', 'M', 'E', 'A', 'N', '1'};
constexpr std::uint32_t kSchemaVersion = 1;

void append_u32(std::vector<std::byte>& bytes, const std::uint32_t value) {
  for (unsigned int shift = 0; shift < 32; shift += 8)
    bytes.push_back(static_cast<std::byte>((value >> shift) & 0xffU));
}

void append_u64(std::vector<std::byte>& bytes, const std::uint64_t value) {
  for (unsigned int shift = 0; shift < 64; shift += 8)
    bytes.push_back(static_cast<std::byte>((value >> shift) & 0xffU));
}

void append_real(std::vector<std::byte>& bytes, const Real value) {
  static_assert(sizeof(Real) == sizeof(std::uint64_t));
  append_u64(bytes, std::bit_cast<std::uint64_t>(value));
}

void atomic_write(const std::filesystem::path& path,
                  const std::span<const std::byte> bytes) {
  auto temporary = path;
  temporary += ".tmp";
  try {
    std::ofstream output(temporary, std::ios::binary | std::ios::trunc);
    if (!output) throw std::runtime_error("unable to open " + temporary.string());
    output.write(reinterpret_cast<const char*>(bytes.data()),
                 static_cast<std::streamsize>(bytes.size()));
    output.close();
    if (!output) throw std::runtime_error("unable to write " + temporary.string());
    std::filesystem::rename(temporary, path);
  } catch (...) {
    std::error_code ignored;
    std::filesystem::remove(temporary, ignored);
    throw;
  }
}

[[nodiscard]] std::string json_escape(const std::string_view value) {
  std::string result;
  for (const char character : value) {
    if (character == '\\' || character == '"') result.push_back('\\');
    result.push_back(character);
  }
  return result;
}

[[nodiscard]] bool selected(const ExperimentConfig& config, const std::string_view id) {
  return config.statistics.fields.empty() ||
         std::ranges::find(config.statistics.fields, id) !=
             config.statistics.fields.end();
}

}  // namespace

std::vector<VisualFieldDescriptor> dry_visual_fields(
    const ExperimentConfig& config, const std::size_t cells, const std::size_t levels,
    const bool has_surface_temperature) {
  std::vector<VisualFieldDescriptor> result;
  std::size_t offset = 0;
  const auto add = [&](const std::string& id, const std::string& unit,
                       const VisualFieldLocation location, const std::size_t count) {
    if (!selected(config, id)) return;
    result.push_back({id, unit, location, offset, count});
    offset += count;
  };
  add("surface_pressure", "Pa", VisualFieldLocation::kSurface, cells);
  add("temperature", "K", VisualFieldLocation::kAtmosphere, cells * levels);
  add("pressure", "Pa", VisualFieldLocation::kAtmosphere, cells * levels);
  add("potential_temperature", "K", VisualFieldLocation::kAtmosphere, cells * levels);
  add("eastward_wind", "m s-1", VisualFieldLocation::kAtmosphere, cells * levels);
  add("northward_wind", "m s-1", VisualFieldLocation::kAtmosphere, cells * levels);
  add("wind_speed", "m s-1", VisualFieldLocation::kAtmosphere, cells * levels);
  if (has_surface_temperature)
    add("surface_temperature", "K", VisualFieldLocation::kSurface, cells);
  for (const auto& tracer : config.tracers)
    add("tracer." + tracer.name, "1", VisualFieldLocation::kAtmosphere, cells * levels);

  if (!config.statistics.fields.empty()) {
    for (const auto& requested : config.statistics.fields) {
      if (std::ranges::none_of(
              result, [&](const auto& field) { return field.id == requested; }))
        throw std::invalid_argument("unsupported statistics field: " + requested);
    }
  }
  if (result.empty()) throw std::invalid_argument("statistics field list is empty");
  return result;
}

std::vector<Real> dry_visual_values(
    const ExperimentConfig& config, const CubedSphereGrid& grid,
    const DryHydrostaticState& state, const DryHydrostaticDerived& derived,
    const std::span<const VisualFieldDescriptor> fields) {
  const auto cells = grid.cell_count();
  const auto levels = derived.levels;
  const auto count = fields.back().value_offset + fields.back().value_count;
  std::vector<Real> result(count);
  const TracerRegistry registry = config.tracers.empty()
                                      ? TracerRegistry::legacy()
                                      : TracerRegistry(config.tracers);
  for (const auto& field : fields) {
    auto output = std::span(result).subspan(field.value_offset, field.value_count);
    if (field.id == "surface_pressure") {
      std::ranges::copy(state.surface_pressure_pa, output.begin());
    } else if (field.id == "temperature") {
      std::ranges::copy(derived.temperature_k, output.begin());
    } else if (field.id == "pressure") {
      std::ranges::copy(derived.pressure_pa, output.begin());
    } else if (field.id == "potential_temperature") {
      std::ranges::copy(derived.potential_temperature_k, output.begin());
    } else if (field.id == "surface_temperature") {
      std::ranges::copy(state.surface_temperature_k, output.begin());
    } else if (field.id.starts_with("tracer.")) {
      const auto name = field.id.substr(7);
      const auto tracer = registry.index_of(name);
      const auto begin = tracer * cells * levels;
      std::ranges::copy(
          std::span(derived.tracer_mixing_ratio).subspan(begin, cells * levels),
          output.begin());
    } else {
      for (std::size_t cell = 0; cell < cells; ++cell) {
        const auto position = grid.cells()[cell].center;
        const Real longitude = std::atan2(position.y, position.x);
        const Vec3 east{-std::sin(longitude), std::cos(longitude), 0.0};
        const Vec3 north = normalize(project_tangent(Vec3{0.0, 0.0, 1.0}, position));
        for (std::size_t level = 0; level < levels; ++level) {
          const auto index = dry_hydrostatic_offset(cell, level, levels);
          const auto velocity = derived.velocity_m_s[index];
          if (field.id == "eastward_wind")
            output[index] = dot(velocity, east);
          else if (field.id == "northward_wind")
            output[index] = dot(velocity, north);
          else if (field.id == "wind_speed")
            output[index] = norm(velocity);
          else
            throw std::logic_error("unknown visual field descriptor");
        }
      }
    }
  }
  for (const auto value : result) require_finite(value, "visual dataset value");
  return result;
}

VisualDatasetWriter::VisualDatasetWriter(std::filesystem::path directory,
                                         const ExperimentConfig& config,
                                         const CubedSphereGrid& grid,
                                         const SurfaceOrography& orography,
                                         const std::span<const Real> land_fraction,
                                         std::vector<VisualFieldDescriptor> fields,
                                         std::string config_fingerprint)
    : directory_(std::move(directory)),
      config_(config),
      grid_(grid),
      orography_(orography),
      land_fraction_(land_fraction.begin(), land_fraction.end()),
      fields_(std::move(fields)),
      config_fingerprint_(std::move(config_fingerprint)) {
  if (!land_fraction_.empty() && land_fraction_.size() != grid_.cell_count())
    throw std::invalid_argument("terrain land fraction shape mismatch");
  std::filesystem::create_directories(directory_ / "means");
  write_terrain();
  write_manifest();
}

std::size_t VisualDatasetWriter::value_count() const noexcept {
  return fields_.empty() ? 0 : fields_.back().value_offset + fields_.back().value_count;
}

void VisualDatasetWriter::write_terrain() const {
  const auto cells = grid_.cell_count();
  const std::uint32_t field_count = land_fraction_.empty() ? 5U : 6U;
  std::vector<std::byte> bytes;
  bytes.reserve(24 + static_cast<std::size_t>(field_count) * cells * sizeof(Real));
  for (const char value : kTerrainMagic) bytes.push_back(static_cast<std::byte>(value));
  append_u32(bytes, kSchemaVersion);
  append_u32(bytes, field_count);
  append_u64(bytes, cells);
  for (const auto& cell : grid_.cells()) append_real(bytes, cell.center.x);
  for (const auto& cell : grid_.cells()) append_real(bytes, cell.center.y);
  for (const auto& cell : grid_.cells()) append_real(bytes, cell.center.z);
  for (const auto value : orography_.surface_geopotential_m2_s2())
    append_real(bytes, value / config_.planet.gravity_m_s2);
  for (const auto value : orography_.surface_geopotential_m2_s2())
    append_real(bytes, value);
  for (const auto value : land_fraction_) append_real(bytes, value);
  atomic_write(directory_ / "terrain.bin", bytes);
}

void VisualDatasetWriter::publish(const diagnostics::PeriodMeanWindow& period) {
  if (period.values.size() != value_count())
    throw std::invalid_argument("period mean does not match visual field layout");
  std::ostringstream filename;
  filename << "period_" << std::setfill('0') << std::setw(6) << period.index << ".bin";
  const auto relative_path = std::string("means/") + filename.str();
  std::vector<std::byte> bytes;
  bytes.reserve(24 + period.values.size() * sizeof(Real));
  for (const char value : kMeanMagic) bytes.push_back(static_cast<std::byte>(value));
  append_u32(bytes, kSchemaVersion);
  append_u32(bytes, static_cast<std::uint32_t>(fields_.size()));
  append_u64(bytes, period.values.size());
  for (const auto value : period.values) append_real(bytes, value);
  atomic_write(directory_ / relative_path, bytes);
  VisualPeriodDescriptor descriptor{.index = period.index,
                                    .scheduled_start_s = period.scheduled_start_s,
                                    .scheduled_end_s = period.scheduled_end_s,
                                    .actual_start_s = period.actual_start_s,
                                    .actual_end_s = period.actual_end_s,
                                    .weight_s = period.weight_s,
                                    .complete = period.complete,
                                    .relative_path = relative_path,
                                    .byte_length = bytes.size()};
  const auto existing = std::ranges::find_if(
      periods_, [&](const auto& value) { return value.index == period.index; });
  if (existing == periods_.end())
    periods_.push_back(std::move(descriptor));
  else
    *existing = std::move(descriptor);
  std::ranges::sort(periods_, {}, &VisualPeriodDescriptor::index);
  write_manifest();
}

void VisualDatasetWriter::set_periods(std::vector<VisualPeriodDescriptor> periods) {
  periods_ = std::move(periods);
  write_manifest();
}

void VisualDatasetWriter::write_manifest() const {
  auto temporary = directory_ / "manifest.json.tmp";
  std::ofstream output(temporary, std::ios::trunc);
  if (!output) throw std::runtime_error("unable to open visual manifest");
  const auto cells = grid_.cell_count();
  const auto levels = static_cast<std::size_t>(config_.vertical.levels);
  output << std::setprecision(std::numeric_limits<Real>::max_digits10)
         << "{\n  \"schemaVersion\": 1,\n  \"modelKind\": \"dry_hydrostatic\",\n"
         << "  \"cellsPerPanel\": " << config_.grid.cells_per_panel << ",\n"
         << "  \"cellCount\": " << cells << ",\n  \"levelCount\": " << levels
         << ",\n  \"arrayOrder\": \"field-cell-level\",\n"
         << "  \"planet\": {\"radiusM\": " << config_.planet.radius_m
         << ", \"gravityMS2\": " << config_.planet.gravity_m_s2 << "},\n"
         << "  \"hybridAHalfPa\": [";
  for (std::size_t i = 0; i < config_.vertical.a_half_pa.size(); ++i) {
    if (i) output << ',';
    output << config_.vertical.a_half_pa[i];
  }
  output << "],\n  \"hybridBHalf\": [";
  for (std::size_t i = 0; i < config_.vertical.b_half.size(); ++i) {
    if (i) output << ',';
    output << config_.vertical.b_half[i];
  }
  output << "],\n  \"configFingerprint\": \"" << json_escape(config_fingerprint_)
         << "\",\n"
         << "  \"terrainSourceFingerprint\": \""
         << json_escape(orography_.source_fingerprint()) << "\",\n"
         << "  \"averagingRule\": \"accepted_step_end_time_weighted_rectangle\",\n"
         << "  \"terrain\": {\"path\": \"terrain.bin\", \"byteLength\": "
         << std::filesystem::file_size(directory_ / "terrain.bin") << ", \"fields\": ["
         << "{\"id\":\"center_x\",\"unit\":\"1\"},"
         << "{\"id\":\"center_y\",\"unit\":\"1\"},"
         << "{\"id\":\"center_z\",\"unit\":\"1\"},"
         << "{\"id\":\"elevation\",\"unit\":\"m\"},"
         << "{\"id\":\"surface_geopotential\",\"unit\":\"m2 s-2\"}";
  if (!land_fraction_.empty()) output << ",{\"id\":\"land_fraction\",\"unit\":\"1\"}";
  output << "]},\n  \"fields\": [";
  for (std::size_t i = 0; i < fields_.size(); ++i) {
    if (i) output << ',';
    const auto& field = fields_[i];
    output << "{\"id\":\"" << json_escape(field.id) << "\",\"unit\":\""
           << json_escape(field.unit) << "\",\"location\":\""
           << (field.location == VisualFieldLocation::kSurface ? "surface"
                                                               : "atmosphere")
           << "\",\"valueOffset\":" << field.value_offset
           << ",\"valueCount\":" << field.value_count << '}';
  }
  output << "],\n  \"periods\": [";
  for (std::size_t i = 0; i < periods_.size(); ++i) {
    if (i) output << ',';
    const auto& period = periods_[i];
    output << "{\"index\":" << period.index
           << ",\"scheduledStartS\":" << period.scheduled_start_s
           << ",\"scheduledEndS\":" << period.scheduled_end_s
           << ",\"actualStartS\":" << period.actual_start_s
           << ",\"actualEndS\":" << period.actual_end_s
           << ",\"weightS\":" << period.weight_s << ",\"coverage\":"
           << period.weight_s / (period.scheduled_end_s - period.scheduled_start_s)
           << ",\"complete\":" << (period.complete ? "true" : "false") << ",\"path\":\""
           << period.relative_path << "\",\"byteLength\":" << period.byte_length << '}';
  }
  output << "]\n}\n";
  output.close();
  if (!output) throw std::runtime_error("unable to write visual manifest");
  std::filesystem::rename(temporary, directory_ / "manifest.json");
}

std::filesystem::path statistics_sidecar_path(
    const std::filesystem::path& checkpoint_path) {
  auto result = checkpoint_path;
  result += ".statistics";
  return result;
}

std::string file_fnv1a64(const std::filesystem::path& path) {
  std::ifstream input(path, std::ios::binary);
  if (!input) throw std::runtime_error("unable to hash file " + path.string());
  std::ostringstream bytes;
  bytes << input.rdbuf();
  if (!input && !input.eof())
    throw std::runtime_error("unable to read file " + path.string());
  return fnv1a64_hex(bytes.str());
}

void write_statistics_sidecar_file(const std::filesystem::path& path,
                                   const StatisticsSidecar& sidecar) {
  auto temporary = path;
  temporary += ".tmp";
  std::ofstream output(temporary, std::ios::trunc);
  if (!output) throw std::runtime_error("unable to open statistics sidecar");
  const auto& state = sidecar.accumulator;
  output.imbue(std::locale::classic());
  output << std::setprecision(std::numeric_limits<Real>::max_digits10)
         << "MPS_STATISTICS_SIDECAR 1\n"
         << sidecar.checkpoint_hash << '\n'
         << sidecar.checkpoint_time_s << '\n'
         << sidecar.checkpoint_step << '\n'
         << state.value_count << '\n'
         << state.start_time_s << '\n'
         << state.period_s << '\n'
         << state.last_step_end_s << '\n'
         << state.has_observation << '\n'
         << state.open_index << '\n'
         << state.open_actual_start_s << '\n'
         << state.open_actual_end_s << '\n'
         << state.open_weight_s << '\n'
         << state.integral.size() << '\n';
  for (const auto value : state.integral) output << value << '\n';
  output << sidecar.periods.size() << '\n';
  for (const auto& period : sidecar.periods)
    output << period.index << ' ' << period.scheduled_start_s << ' '
           << period.scheduled_end_s << ' ' << period.actual_start_s << ' '
           << period.actual_end_s << ' ' << period.weight_s << ' ' << period.complete
           << ' ' << period.byte_length << ' ' << std::quoted(period.relative_path)
           << '\n';
  output.close();
  if (!output) throw std::runtime_error("unable to write statistics sidecar");
  std::filesystem::rename(temporary, path);
}

StatisticsSidecar read_statistics_sidecar_file(
    const std::filesystem::path& path, const std::string_view expected_checkpoint_hash,
    const Real expected_time_s, const std::uint64_t expected_step,
    const std::size_t expected_value_count, const Real expected_start_time_s,
    const Real expected_period_s) {
  std::ifstream input(path);
  if (!input) throw std::runtime_error("unable to open statistics sidecar");
  input.imbue(std::locale::classic());
  std::string magic;
  unsigned int version = 0;
  StatisticsSidecar result;
  auto& state = result.accumulator;
  std::size_t integral_size = 0;
  input >> magic >> version >> result.checkpoint_hash >> result.checkpoint_time_s >>
      result.checkpoint_step >> state.value_count >> state.start_time_s >>
      state.period_s >> state.last_step_end_s >> state.has_observation >>
      state.open_index >> state.open_actual_start_s >> state.open_actual_end_s >>
      state.open_weight_s >> integral_size;
  if (!input || magic != "MPS_STATISTICS_SIDECAR" || version != 1)
    throw std::runtime_error("invalid statistics sidecar header");
  if (integral_size > 100'000'000 || integral_size != state.value_count)
    throw std::runtime_error("invalid statistics sidecar value count");
  state.integral.resize(integral_size);
  for (auto& value : state.integral) input >> value;
  std::size_t period_count = 0;
  input >> period_count;
  if (!input || period_count > 1'000'000)
    throw std::runtime_error("invalid statistics sidecar period count");
  result.periods.resize(period_count);
  for (auto& period : result.periods)
    input >> period.index >> period.scheduled_start_s >> period.scheduled_end_s >>
        period.actual_start_s >> period.actual_end_s >> period.weight_s >>
        period.complete >> period.byte_length >> std::quoted(period.relative_path);
  input >> std::ws;
  if (!input.eof()) throw std::runtime_error("statistics sidecar has trailing data");
  if (result.checkpoint_hash != expected_checkpoint_hash ||
      result.checkpoint_time_s != expected_time_s ||
      result.checkpoint_step != expected_step)
    throw std::runtime_error("statistics sidecar checkpoint mismatch");
  if (state.value_count != expected_value_count ||
      state.start_time_s != expected_start_time_s ||
      state.period_s != expected_period_s)
    throw std::runtime_error("statistics sidecar configuration mismatch");
  // Reuse accumulator validation for all numeric state.
  static_cast<void>(diagnostics::PeriodMeanAccumulator(state));
  return result;
}

}  // namespace mps
