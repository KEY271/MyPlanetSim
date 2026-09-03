#include "myplanetsim/geometry/cubed_sphere.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

namespace mps {

std::string_view panel_name(const Panel panel) noexcept {
  constexpr std::array<std::string_view, 6> names{"PX", "PY", "NX", "NY", "PZ", "NZ"};
  return names[panel_index(panel)];
}

PanelBasis panel_basis(const Panel panel) noexcept {
  constexpr std::array<PanelBasis, 6> bases{{
      {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}},
      {{0, 1, 0}, {-1, 0, 0}, {0, 0, 1}},
      {{-1, 0, 0}, {0, -1, 0}, {0, 0, 1}},
      {{0, -1, 0}, {1, 0, 0}, {0, 0, 1}},
      {{0, 0, 1}, {1, 0, 0}, {0, 1, 0}},
      {{0, 0, -1}, {1, 0, 0}, {0, -1, 0}},
  }};
  return bases[panel_index(panel)];
}

Vec3 map_to_unit_sphere(const Panel panel, const Real alpha, const Real beta) {
  if (!std::isfinite(alpha) || !std::isfinite(beta)) {
    throw std::invalid_argument("panel coordinates must be finite");
  }
  const auto basis = panel_basis(panel);
  return normalize(basis.normal + std::tan(alpha) * basis.alpha +
                   std::tan(beta) * basis.beta);
}

PanelCoordinates inverse_map(const Vec3 position) {
  const Vec3 point = normalize(position);
  const std::array<Real, 3> magnitudes{std::abs(point.x), std::abs(point.y),
                                       std::abs(point.z)};
  const std::size_t axis = static_cast<std::size_t>(std::distance(
      magnitudes.begin(), std::max_element(magnitudes.begin(), magnitudes.end())));
  Panel panel = Panel::kPositiveX;
  if (axis == 0) {
    panel = point.x >= 0.0 ? Panel::kPositiveX : Panel::kNegativeX;
  } else if (axis == 1) {
    panel = point.y >= 0.0 ? Panel::kPositiveY : Panel::kNegativeY;
  } else {
    panel = point.z >= 0.0 ? Panel::kPositiveZ : Panel::kNegativeZ;
  }
  const auto basis = panel_basis(panel);
  const Real denominator = dot(point, basis.normal);
  return {panel, std::atan2(dot(point, basis.alpha), denominator),
          std::atan2(dot(point, basis.beta), denominator)};
}

TangentBasis tangent_basis(const Panel panel, const Real alpha, const Real beta) {
  const auto fixed = panel_basis(panel);
  const Vec3 position = map_to_unit_sphere(panel, alpha, beta);
  const Vec3 alpha_direction = normalize(project_tangent(fixed.alpha, position));
  Vec3 beta_direction = normalize(cross(position, alpha_direction));
  if (dot(beta_direction, fixed.beta) < 0.0) {
    beta_direction = -beta_direction;
  }
  return {alpha_direction, beta_direction};
}

Vec3 from_tangent_components(const TangentComponents components,
                             const TangentBasis& basis) noexcept {
  return components.alpha * basis.alpha + components.beta * basis.beta;
}

TangentComponents to_tangent_components(const Vec3 vector, const TangentBasis& basis) {
  if (!is_finite(vector)) {
    throw std::invalid_argument("tangent vector must be finite");
  }
  return {dot(vector, basis.alpha), dot(vector, basis.beta)};
}

}  // namespace mps
