#include "myplanetsim/io/frame.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstring>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <vector>

#include "myplanetsim/core/validation.hpp"

namespace mps {
namespace {

constexpr std::array<unsigned char, 8> kMagic{'M', 'P', 'S', 'F', 'R', 'A', 'M', '1'};
constexpr std::array<unsigned char, 8> kMagicV2{'M', 'P', 'S', 'F', 'R', 'A', 'M', '2'};
constexpr std::uint32_t kVersion = 1;
constexpr std::uint64_t kMaximumPayloadBytes = 1024ULL * 1024ULL * 1024ULL;

void append_u32(std::vector<unsigned char>& bytes, const std::uint32_t value) {
  for (unsigned int shift = 0; shift < 32; shift += 8) {
    bytes.push_back(static_cast<unsigned char>(value >> shift));
  }
}

void append_u64(std::vector<unsigned char>& bytes, const std::uint64_t value) {
  for (unsigned int shift = 0; shift < 64; shift += 8) {
    bytes.push_back(static_cast<unsigned char>(value >> shift));
  }
}

[[nodiscard]] std::uint32_t read_u32(const std::vector<unsigned char>& bytes,
                                     std::size_t& offset) {
  if (offset + 4 > bytes.size()) {
    throw std::runtime_error("frame is truncated");
  }
  std::uint32_t value = 0;
  for (unsigned int shift = 0; shift < 32; shift += 8) {
    value |= static_cast<std::uint32_t>(bytes[offset++]) << shift;
  }
  return value;
}

[[nodiscard]] std::uint64_t read_u64(const std::vector<unsigned char>& bytes,
                                     std::size_t& offset) {
  if (offset + 8 > bytes.size()) {
    throw std::runtime_error("frame is truncated");
  }
  std::uint64_t value = 0;
  for (unsigned int shift = 0; shift < 64; shift += 8) {
    value |= static_cast<std::uint64_t>(bytes[offset++]) << shift;
  }
  return value;
}

void append_real(std::vector<unsigned char>& bytes, const Real value) {
  append_u64(bytes, std::bit_cast<std::uint64_t>(value));
}

[[nodiscard]] Real read_real(const std::vector<unsigned char>& bytes,
                             std::size_t& offset) {
  return std::bit_cast<Real>(read_u64(bytes, offset));
}

[[nodiscard]] std::vector<unsigned char> read_bytes(const std::filesystem::path& path) {
  std::ifstream input(path, std::ios::binary | std::ios::ate);
  if (!input) {
    throw std::runtime_error("unable to open frame: " + path.string());
  }
  const auto size = input.tellg();
  if (size < 0 || static_cast<std::uint64_t>(size) > kMaximumPayloadBytes) {
    throw std::runtime_error("frame is outside the supported size");
  }
  std::vector<unsigned char> bytes(static_cast<std::size_t>(size));
  input.seekg(0);
  if (!bytes.empty()) {
    input.read(reinterpret_cast<char*>(bytes.data()),
               static_cast<std::streamsize>(bytes.size()));
  }
  if (!input) {
    throw std::runtime_error("failed while reading frame");
  }
  return bytes;
}

}  // namespace

void write_frame_file(const std::filesystem::path& path,
                      const FrameV1Metadata& metadata, const ShallowWaterState& state) {
  if (metadata.schema_version != kVersion || metadata.cells_per_panel <= 0 ||
      metadata.cells_per_panel > 1024 || metadata.config_fingerprint.empty() ||
      metadata.config_fingerprint.size() > 128) {
    throw std::invalid_argument("invalid FrameV1 metadata");
  }
  require_non_negative(metadata.time_s, "frame time_s");
  if (state.depth.empty() || state.momentum.size() != state.depth.size()) {
    throw std::invalid_argument("frame state shape is invalid");
  }
  const auto n = static_cast<std::uint64_t>(metadata.cells_per_panel);
  const auto expected_cells = 6 * n * n;
  if (state.depth.size() != expected_cells) {
    throw std::invalid_argument("frame state size does not match cells_per_panel");
  }
  std::vector<unsigned char> bytes;
  bytes.reserve(64 + metadata.config_fingerprint.size() + 32 * state.depth.size());
  bytes.insert(bytes.end(), kMagic.begin(), kMagic.end());
  append_u32(bytes, kVersion);
  append_u32(bytes, 0);
  append_u64(bytes, n);
  append_real(bytes, metadata.time_s);
  append_u64(bytes, metadata.step);
  append_u64(bytes, state.depth.size());
  append_u32(bytes, static_cast<std::uint32_t>(metadata.config_fingerprint.size()));
  bytes.insert(bytes.end(), metadata.config_fingerprint.begin(),
               metadata.config_fingerprint.end());
  const auto flat = flatten_shallow_water_state(state);
  for (const Real value : flat) {
    require_finite(value, "frame value");
    append_real(bytes, value);
  }

  if (!path.parent_path().empty()) {
    std::filesystem::create_directories(path.parent_path());
  }
  auto temporary_path = path;
  temporary_path += ".tmp";
  std::ofstream output(temporary_path, std::ios::binary | std::ios::trunc);
  if (!output) {
    throw std::runtime_error("unable to open temporary frame: " +
                             temporary_path.string());
  }
  output.write(reinterpret_cast<const char*>(bytes.data()),
               static_cast<std::streamsize>(bytes.size()));
  output.close();
  if (!output) {
    throw std::runtime_error("failed while writing frame");
  }
  std::filesystem::rename(temporary_path, path);
}

FrameV1 read_frame_file(const std::filesystem::path& path,
                        const std::string_view expected_config_fingerprint) {
  const auto bytes = read_bytes(path);
  std::size_t offset = 0;
  if (bytes.size() < kMagic.size() ||
      !std::equal(kMagic.begin(), kMagic.end(), bytes.begin())) {
    throw std::runtime_error("invalid FrameV1 magic");
  }
  offset = kMagic.size();
  if (read_u32(bytes, offset) != kVersion || read_u32(bytes, offset) != 0) {
    throw std::runtime_error("unsupported or invalid FrameV1 header");
  }
  const auto cells_per_panel = read_u64(bytes, offset);
  if (cells_per_panel == 0 || cells_per_panel > 1024) {
    throw std::runtime_error("FrameV1 cells_per_panel is outside the supported range");
  }
  const Real time_s = read_real(bytes, offset);
  const auto step = read_u64(bytes, offset);
  const auto cell_count = read_u64(bytes, offset);
  const auto fingerprint_size = read_u32(bytes, offset);
  if (fingerprint_size == 0 || fingerprint_size > 128 ||
      offset + fingerprint_size > bytes.size()) {
    throw std::runtime_error("invalid FrameV1 fingerprint");
  }
  const std::string fingerprint(reinterpret_cast<const char*>(bytes.data() + offset),
                                fingerprint_size);
  offset += fingerprint_size;
  if (!expected_config_fingerprint.empty() &&
      fingerprint != expected_config_fingerprint) {
    throw std::runtime_error("frame configuration fingerprint mismatch");
  }
  const auto expected_cells = 6 * cells_per_panel * cells_per_panel;
  if (cell_count != expected_cells ||
      cell_count > std::numeric_limits<std::size_t>::max() / 4) {
    throw std::runtime_error("FrameV1 cell count does not match shape");
  }
  const auto value_count = static_cast<std::size_t>(4 * cell_count);
  if (offset + value_count * sizeof(Real) != bytes.size()) {
    throw std::runtime_error("FrameV1 payload has the wrong size");
  }
  std::vector<Real> flat;
  flat.reserve(value_count);
  for (std::size_t index = 0; index < value_count; ++index) {
    const Real value = read_real(bytes, offset);
    require_finite(value, "frame value");
    flat.push_back(value);
  }
  return {.metadata = {.schema_version = kVersion,
                       .cells_per_panel = static_cast<Index>(cells_per_panel),
                       .time_s = time_s,
                       .step = step,
                       .config_fingerprint = fingerprint},
          .state = unflatten_shallow_water_state(time_s, step, flat,
                                                 static_cast<std::size_t>(cell_count)),
          .byte_length = static_cast<std::uint64_t>(bytes.size())};
}

void write_frame_v2_file(const std::filesystem::path& path,const FrameV2Metadata& m,const std::span<const Real> ps,const DryHydrostaticDerived& d){const auto n=static_cast<std::uint64_t>(m.cells_per_panel),cells=6*n*n,levels=static_cast<std::uint64_t>(m.levels),volume=cells*levels;if(m.schema_version!=2||n==0||n>24||levels==0||levels>30||ps.size()!=cells||d.cells!=cells||d.levels!=levels||d.pressure_pa.size()!=volume||d.potential_temperature_k.size()!=volume||d.temperature_k.size()!=volume||d.velocity_m_s.size()!=volume||d.tracer_mixing_ratio.size()!=volume||m.config_fingerprint.empty()||m.config_fingerprint.size()>128)throw std::invalid_argument("invalid FrameV2 shape or metadata");std::vector<unsigned char>b;b.insert(b.end(),kMagicV2.begin(),kMagicV2.end());append_u32(b,2);append_u32(b,0);append_u64(b,n);append_u64(b,levels);append_real(b,m.time_s);append_u64(b,m.step);append_u64(b,cells);append_u32(b,static_cast<std::uint32_t>(m.config_fingerprint.size()));b.insert(b.end(),m.config_fingerprint.begin(),m.config_fingerprint.end());auto add=[&](Real v){require_finite(v,"FrameV2 value");append_real(b,v);};for(auto v:ps)add(v);for(auto v:d.pressure_pa)add(v);for(auto v:d.potential_temperature_k)add(v);for(auto v:d.temperature_k)add(v);for(auto v:d.velocity_m_s)add(v.x);for(auto v:d.velocity_m_s)add(v.y);for(auto v:d.velocity_m_s)add(v.z);for(auto v:d.tracer_mixing_ratio)add(v);if(b.size()>256ULL*1024ULL*1024ULL)throw std::runtime_error("FrameV2 exceeds byte budget");if(!path.parent_path().empty())std::filesystem::create_directories(path.parent_path());auto tmp=path;tmp+=".tmp";std::ofstream out(tmp,std::ios::binary|std::ios::trunc);out.write(reinterpret_cast<const char*>(b.data()),static_cast<std::streamsize>(b.size()));out.close();if(!out)throw std::runtime_error("failed while writing FrameV2");std::filesystem::rename(tmp,path);}

FrameV2 read_frame_v2_file(const std::filesystem::path& path,const std::string_view expected){auto b=read_bytes(path);std::size_t o=0;if(b.size()<kMagicV2.size()||!std::equal(kMagicV2.begin(),kMagicV2.end(),b.begin()))throw std::runtime_error("invalid FrameV2 magic");o=8;if(read_u32(b,o)!=2||read_u32(b,o)!=0)throw std::runtime_error("invalid FrameV2 header");auto n=read_u64(b,o),K=read_u64(b,o);auto time=read_real(b,o);auto step=read_u64(b,o),C=read_u64(b,o);auto fs=read_u32(b,o);if(n==0||n>24||K==0||K>30||C!=6*n*n||fs==0||fs>128||o+fs>b.size())throw std::runtime_error("invalid FrameV2 shape");std::string fingerprint(reinterpret_cast<const char*>(b.data()+o),fs);o+=fs;if(!expected.empty()&&fingerprint!=expected)throw std::runtime_error("frame configuration fingerprint mismatch");auto volume=C*K,value_count=C+7*volume;if(o+value_count*sizeof(Real)!=b.size())throw std::runtime_error("FrameV2 payload has wrong size");auto read_field=[&](std::size_t count){std::vector<Real> values;values.reserve(count);for(std::size_t i=0;i<count;++i){auto v=read_real(b,o);require_finite(v,"FrameV2 value");values.push_back(v);}return values;};FrameV2 f;f.metadata={2,static_cast<Index>(n),static_cast<Index>(K),time,step,fingerprint};f.surface_pressure_pa=read_field(C);f.derived.cells=C;f.derived.levels=K;f.derived.pressure_pa=read_field(volume);f.derived.potential_temperature_k=read_field(volume);f.derived.temperature_k=read_field(volume);auto x=read_field(volume),y=read_field(volume),z=read_field(volume);f.derived.velocity_m_s.resize(volume);for(std::size_t i=0;i<volume;++i)f.derived.velocity_m_s[i]={x[i],y[i],z[i]};f.derived.tracer_mixing_ratio=read_field(volume);f.byte_length=b.size();return f;}

}  // namespace mps
