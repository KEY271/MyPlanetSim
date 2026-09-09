#pragma once

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#include "myplanetsim/config/experiment_config.hpp"
#include "myplanetsim/diagnostics/period_mean.hpp"
#include "myplanetsim/dynamics/dry_hydrostatic_state.hpp"
#include "myplanetsim/dynamics/surface_orography.hpp"
#include "myplanetsim/grid/cubed_sphere_grid.hpp"

namespace mps {

enum class VisualFieldLocation { kSurface, kAtmosphere };

struct VisualFieldDescriptor {
  std::string id;
  std::string unit;
  VisualFieldLocation location = VisualFieldLocation::kSurface;
  std::size_t value_offset = 0;
  std::size_t value_count = 0;
};

struct VisualPeriodDescriptor {
  std::uint64_t index = 0;
  Real scheduled_start_s = 0.0;
  Real scheduled_end_s = 0.0;
  Real actual_start_s = 0.0;
  Real actual_end_s = 0.0;
  Real weight_s = 0.0;
  bool complete = false;
  std::string relative_path;
  std::uintmax_t byte_length = 0;
};

[[nodiscard]] std::vector<VisualFieldDescriptor> dry_visual_fields(
    const ExperimentConfig& config, std::size_t cells, std::size_t levels,
    bool has_surface_temperature);
[[nodiscard]] std::vector<Real> dry_visual_values(
    const ExperimentConfig& config, const CubedSphereGrid& grid,
    const DryHydrostaticState& state, const DryHydrostaticDerived& derived,
    std::span<const VisualFieldDescriptor> fields);

class VisualDatasetWriter {
 public:
  VisualDatasetWriter(std::filesystem::path directory, const ExperimentConfig& config,
                      const CubedSphereGrid& grid, const SurfaceOrography& orography,
                      std::span<const Real> land_fraction,
                      std::vector<VisualFieldDescriptor> fields,
                      std::string config_fingerprint);

  void publish(const diagnostics::PeriodMeanWindow& period);
  void set_periods(std::vector<VisualPeriodDescriptor> periods);
  [[nodiscard]] const std::vector<VisualPeriodDescriptor>& periods() const noexcept {
    return periods_;
  }
  [[nodiscard]] const std::vector<VisualFieldDescriptor>& fields() const noexcept {
    return fields_;
  }
  [[nodiscard]] std::size_t value_count() const noexcept;

 private:
  std::filesystem::path directory_;
  const ExperimentConfig& config_;
  const CubedSphereGrid& grid_;
  const SurfaceOrography& orography_;
  std::vector<Real> land_fraction_;
  std::vector<VisualFieldDescriptor> fields_;
  std::string config_fingerprint_;
  std::vector<VisualPeriodDescriptor> periods_;

  void write_terrain() const;
  void write_manifest() const;
};

struct StatisticsSidecar {
  std::string checkpoint_hash;
  Real checkpoint_time_s = 0.0;
  std::uint64_t checkpoint_step = 0;
  diagnostics::PeriodMeanAccumulatorState accumulator;
  std::vector<VisualPeriodDescriptor> periods;
};

[[nodiscard]] std::filesystem::path statistics_sidecar_path(
    const std::filesystem::path& checkpoint_path);
[[nodiscard]] std::string file_fnv1a64(const std::filesystem::path& path);
void write_statistics_sidecar_file(const std::filesystem::path& path,
                                   const StatisticsSidecar& sidecar);
[[nodiscard]] StatisticsSidecar read_statistics_sidecar_file(
    const std::filesystem::path& path, std::string_view expected_checkpoint_hash,
    Real expected_time_s, std::uint64_t expected_step, std::size_t expected_value_count,
    Real expected_start_time_s, Real expected_period_s);

}  // namespace mps
