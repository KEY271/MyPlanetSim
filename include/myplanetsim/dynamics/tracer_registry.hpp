#pragma once

#include <cstddef>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#include "myplanetsim/core/types.hpp"

namespace mps {

enum class TracerRole { kPassive, kWaterVapor };

struct TracerDescriptor {
  std::string name{};
  TracerRole role = TracerRole::kPassive;
  Real initial_mixing_ratio = 0.0;
  std::optional<Real> initial_relative_humidity{};
  bool require_nonnegative = true;
  bool horizontal_diffusion = true;
};

class TracerRegistry {
 public:
  explicit TracerRegistry(std::vector<TracerDescriptor> tracers);

  [[nodiscard]] static TracerRegistry legacy();
  [[nodiscard]] std::size_t size() const noexcept { return tracers_.size(); }
  [[nodiscard]] bool empty() const noexcept { return tracers_.empty(); }
  [[nodiscard]] const TracerDescriptor& operator[](std::size_t index) const;
  [[nodiscard]] std::span<const TracerDescriptor> tracers() const noexcept {
    return tracers_;
  }
  [[nodiscard]] std::size_t index_of(std::string_view name) const;
  [[nodiscard]] std::optional<std::size_t> water_vapor_index() const noexcept;

 private:
  std::vector<TracerDescriptor> tracers_;
  std::optional<std::size_t> water_vapor_index_;
};

[[nodiscard]] std::string_view tracer_role_name(TracerRole role) noexcept;

}  // namespace mps
