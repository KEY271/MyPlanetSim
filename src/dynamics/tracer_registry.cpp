#include "myplanetsim/dynamics/tracer_registry.hpp"

#include <algorithm>
#include <cctype>
#include <stdexcept>

namespace mps {
namespace {

[[nodiscard]] bool valid_name(const std::string_view name) {
  if (name.empty() ||
      !(std::isalpha(static_cast<unsigned char>(name.front())) || name.front() == '_'))
    return false;
  return std::ranges::all_of(name, [](const char value) {
    const auto byte = static_cast<unsigned char>(value);
    return std::isalnum(byte) || value == '_';
  });
}

}  // namespace

TracerRegistry::TracerRegistry(std::vector<TracerDescriptor> tracers)
    : tracers_(std::move(tracers)) {
  if (tracers_.empty())
    throw std::invalid_argument("tracer registry must not be empty");
  for (std::size_t index = 0; index < tracers_.size(); ++index) {
    const auto& tracer = tracers_[index];
    if (!valid_name(tracer.name))
      throw std::invalid_argument("tracer name must be an identifier");
    for (std::size_t previous = 0; previous < index; ++previous)
      if (tracers_[previous].name == tracer.name)
        throw std::invalid_argument("tracer names must be unique");
    if (tracer.role == TracerRole::kWaterVapor) {
      if (water_vapor_index_.has_value())
        throw std::invalid_argument("at most one water-vapor tracer is allowed");
      water_vapor_index_ = index;
    }
  }
}

TracerRegistry TracerRegistry::legacy() {
  return TracerRegistry({TracerDescriptor{.name = "passive"}});
}

const TracerDescriptor& TracerRegistry::operator[](const std::size_t index) const {
  if (index >= tracers_.size()) throw std::out_of_range("tracer index is out of range");
  return tracers_[index];
}

std::size_t TracerRegistry::index_of(const std::string_view name) const {
  for (std::size_t index = 0; index < tracers_.size(); ++index)
    if (tracers_[index].name == name) return index;
  throw std::out_of_range("unknown tracer name");
}

std::optional<std::size_t> TracerRegistry::water_vapor_index() const noexcept {
  return water_vapor_index_;
}

std::string_view tracer_role_name(const TracerRole role) noexcept {
  switch (role) {
    case TracerRole::kPassive:
      return "passive";
    case TracerRole::kWaterVapor:
      return "water_vapor";
  }
  return "unknown";
}

}  // namespace mps
