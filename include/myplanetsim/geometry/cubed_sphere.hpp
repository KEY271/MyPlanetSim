#pragma once

#include <array>
#include <string_view>

#include "myplanetsim/geometry/vec3.hpp"

namespace mps {

enum class Panel {
  kPositiveX,
  kPositiveY,
  kNegativeX,
  kNegativeY,
  kPositiveZ,
  kNegativeZ
};

inline constexpr std::array<Panel, 6> kPanels{Panel::kPositiveX, Panel::kPositiveY,
                                              Panel::kNegativeX, Panel::kNegativeY,
                                              Panel::kPositiveZ, Panel::kNegativeZ};

struct PanelBasis {
  Vec3 normal;
  Vec3 alpha;
  Vec3 beta;
};

struct PanelCoordinates {
  Panel panel;
  Real alpha;
  Real beta;
};

struct TangentBasis {
  Vec3 alpha;
  Vec3 beta;
};

struct TangentComponents {
  Real alpha;
  Real beta;
};

[[nodiscard]] constexpr std::size_t panel_index(Panel panel) noexcept {
  return static_cast<std::size_t>(panel);
}
[[nodiscard]] std::string_view panel_name(Panel panel) noexcept;
[[nodiscard]] PanelBasis panel_basis(Panel panel) noexcept;
[[nodiscard]] Vec3 map_to_unit_sphere(Panel panel, Real alpha, Real beta);
[[nodiscard]] PanelCoordinates inverse_map(Vec3 position);
[[nodiscard]] TangentBasis tangent_basis(Panel panel, Real alpha, Real beta);
[[nodiscard]] Vec3 from_tangent_components(TangentComponents components,
                                           const TangentBasis& basis) noexcept;
[[nodiscard]] TangentComponents to_tangent_components(Vec3 vector,
                                                      const TangentBasis& basis);

}  // namespace mps
