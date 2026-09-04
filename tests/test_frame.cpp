#include <bit>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <vector>

#include "myplanetsim/io/frame.hpp"
#include "support/test.hpp"

namespace {

[[nodiscard]] mps::ShallowWaterState state(const mps::CubedSphereGrid& grid) {
  mps::ShallowWaterState value{.time_s = 3.5,
                               .step = 7,
                               .depth = std::vector<mps::Real>(grid.cell_count()),
                               .momentum = std::vector<mps::Vec3>(grid.cell_count())};
  for (std::size_t index = 0; index < grid.cell_count(); ++index) {
    value.depth[index] = 2.0 + static_cast<mps::Real>(index);
    value.momentum[index] =
        mps::project_tangent({1.0, -2.0, 3.0}, grid.cells()[index].center);
  }
  return value;
}

[[nodiscard]] std::filesystem::path path() {
  return std::filesystem::temp_directory_path() / "myplanetsim_frame_test.bin";
}

}  // namespace

MPS_TEST_CASE("FrameV1 round trips binary state bitwise") {
  const mps::CubedSphereGrid grid(2, 1.0);
  const auto original = state(grid);
  const auto output = path();
  mps::write_frame_file(output,
                        {.cells_per_panel = 2,
                         .time_s = original.time_s,
                         .step = original.step,
                         .config_fingerprint = "0123456789abcdef"},
                        original);
  const auto restored = mps::read_frame_file(output, "0123456789abcdef");
  MPS_CHECK_EQ(restored.metadata.cells_per_panel, 2);
  MPS_CHECK_EQ(restored.state.step, original.step);
  MPS_CHECK_EQ(restored.byte_length,
               static_cast<std::uint64_t>(std::filesystem::file_size(output)));
  for (std::size_t index = 0; index < original.depth.size(); ++index) {
    MPS_CHECK_EQ(std::bit_cast<std::uint64_t>(restored.state.depth[index]),
                 std::bit_cast<std::uint64_t>(original.depth[index]));
  }
  std::filesystem::remove(output);
}

MPS_TEST_CASE("FrameV1 rejects truncation and wrong fingerprint") {
  const mps::CubedSphereGrid grid(1, 1.0);
  const auto output = path();
  mps::write_frame_file(output,
                        {.cells_per_panel = 1,
                         .time_s = 0.0,
                         .step = 0,
                         .config_fingerprint = "fingerprint"},
                        state(grid));
  MPS_CHECK_THROWS_AS(mps::read_frame_file(output, "other"), std::runtime_error);
  const auto truncated = output.string() + ".truncated";
  {
    std::ifstream input(output, std::ios::binary);
    std::ofstream copy(truncated, std::ios::binary | std::ios::trunc);
    std::vector<char> bytes((std::istreambuf_iterator<char>(input)), {});
    copy.write(bytes.data(), static_cast<std::streamsize>(bytes.size() - 1));
  }
  MPS_CHECK_THROWS_AS(mps::read_frame_file(truncated), std::runtime_error);
  std::filesystem::remove(output);
  std::filesystem::remove(truncated);
}

MPS_TEST_CASE("FrameV2 round trips fixed cell-major level fields") {const mps::CubedSphereGrid grid(1,1);mps::DryHydrostaticDerived d{.cells=6,.levels=2,.pressure_pa=std::vector<mps::Real>(12,50000),.velocity_m_s=std::vector<mps::Vec3>(12,{1,2,3}),.potential_temperature_k=std::vector<mps::Real>(12,300),.tracer_mixing_ratio=std::vector<mps::Real>(12,.5),.temperature_k=std::vector<mps::Real>(12,280)};auto output=path();mps::write_frame_v2_file(output,{.cells_per_panel=1,.levels=2,.time_s=3,.step=4,.config_fingerprint="fp"},std::vector<mps::Real>(6,100000),d);auto f=mps::read_frame_v2_file(output,"fp");MPS_CHECK_EQ(f.metadata.levels,2);MPS_CHECK_EQ(f.derived.pressure_pa.size(),12U);MPS_CHECK_EQ(f.derived.velocity_m_s[3].z,3.0);std::filesystem::remove(output);}

int main() { return mps::test::run_all(); }
