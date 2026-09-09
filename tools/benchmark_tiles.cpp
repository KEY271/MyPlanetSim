#include <bit>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string_view>

#include "myplanetsim/config/experiment_config.hpp"
#include "myplanetsim/dynamics/dry_hydrostatic_driver.hpp"

namespace {

[[nodiscard]] std::size_t parse_count(const std::string_view text) {
  std::size_t result = 0;
  for (const char digit : text) {
    if (digit < '0' || digit > '9') throw std::invalid_argument("invalid RHS count");
    result = 10 * result + static_cast<std::size_t>(digit - '0');
  }
  if (result == 0) throw std::invalid_argument("RHS count must be positive");
  return result;
}

void mix(const mps::Real value, std::uint64_t& hash) {
  hash ^= std::bit_cast<std::uint64_t>(value);
  hash *= 1099511628211ULL;
}

[[nodiscard]] std::uint64_t rhs_hash(const mps::DryHydrostaticRhs& rhs) {
  std::uint64_t hash = 1469598103934665603ULL;
  for (const auto value : rhs.surface_pressure_pa_s) mix(value, hash);
  for (const auto value : rhs.tendency.air_mass) mix(value, hash);
  for (const auto value : rhs.tendency.momentum) {
    mix(value.x, hash);
    mix(value.y, hash);
    mix(value.z, hash);
  }
  for (const auto value : rhs.tendency.potential_temperature_mass) mix(value, hash);
  for (const auto value : rhs.tendency.tracer_mass) mix(value, hash);
  return hash;
}

}  // namespace

int main(const int argc, char** argv) {
  try {
    if (argc != 2 && argc != 4)
      throw std::invalid_argument("usage: benchmark_tiles CONFIG [--rhs-count COUNT]");
    const auto rhs_count = argc == 4 && std::string_view(argv[2]) == "--rhs-count"
                               ? parse_count(argv[3])
                               : 3U;
    if (argc == 4 && std::string_view(argv[2]) != "--rhs-count")
      throw std::invalid_argument("usage: benchmark_tiles CONFIG [--rhs-count COUNT]");
    const auto config = mps::load_experiment_config(argv[1]);
    std::uint64_t reference_hash = 0;
    for (const std::size_t width : {8U, 16U, 32U}) {
      const mps::DryHydrostaticDriver driver(config, width);
      driver.enable_rhs_profiling(true);
      const auto state = driver.initial_state();
      mps::DryHydrostaticRhs rhs;
      const auto start = std::chrono::steady_clock::now();
      for (std::size_t repeat = 0; repeat < rhs_count; ++repeat) driver.rhs(state, rhs);
      const auto elapsed =
          std::chrono::duration<double>(std::chrono::steady_clock::now() - start)
              .count();
      const auto hash = rhs_hash(rhs);
      if (width == 8)
        reference_hash = hash;
      else if (hash != reference_hash)
        throw std::runtime_error("tile width changed the RHS bit pattern");
      const auto& profile = driver.rhs_profile();
      std::cout << "tile_width=" << width
                << " seconds_per_rhs=" << elapsed / static_cast<double>(rhs_count)
                << " tile_count=" << profile.tile_count
                << " halo_cell_references=" << profile.tile_halo_cell_references
                << " metadata_bytes=" << profile.tile_metadata_bytes
                << " rhs_hash=" << hash << '\n';
    }
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 2;
  }
}
