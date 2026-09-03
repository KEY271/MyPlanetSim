#pragma once

#include <cmath>
#include <cstddef>
#include <stdexcept>

#include "myplanetsim/core/types.hpp"

namespace mps {

struct Vec3 {
  Real x = 0.0;
  Real y = 0.0;
  Real z = 0.0;

  [[nodiscard]] constexpr Real& operator[](const std::size_t index) {
    if (index > 2) {
      throw std::out_of_range("Vec3 component index is out of range");
    }
    return index == 0 ? x : (index == 1 ? y : z);
  }

  [[nodiscard]] constexpr const Real& operator[](const std::size_t index) const {
    if (index > 2) {
      throw std::out_of_range("Vec3 component index is out of range");
    }
    return index == 0 ? x : (index == 1 ? y : z);
  }
};

[[nodiscard]] constexpr Vec3 operator+(const Vec3 a, const Vec3 b) noexcept {
  return {a.x + b.x, a.y + b.y, a.z + b.z};
}
[[nodiscard]] constexpr Vec3 operator-(const Vec3 a, const Vec3 b) noexcept {
  return {a.x - b.x, a.y - b.y, a.z - b.z};
}
[[nodiscard]] constexpr Vec3 operator-(const Vec3 a) noexcept {
  return {-a.x, -a.y, -a.z};
}
[[nodiscard]] constexpr Vec3 operator*(const Vec3 a, const Real scale) noexcept {
  return {a.x * scale, a.y * scale, a.z * scale};
}
[[nodiscard]] constexpr Vec3 operator*(const Real scale, const Vec3 a) noexcept {
  return a * scale;
}
[[nodiscard]] constexpr Vec3 operator/(const Vec3 a, const Real scale) noexcept {
  return {a.x / scale, a.y / scale, a.z / scale};
}
[[nodiscard]] constexpr Real dot(const Vec3 a, const Vec3 b) noexcept {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}
[[nodiscard]] constexpr Vec3 cross(const Vec3 a, const Vec3 b) noexcept {
  return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
[[nodiscard]] constexpr Real norm_squared(const Vec3 value) noexcept {
  return dot(value, value);
}
[[nodiscard]] inline Real norm(const Vec3 value) noexcept {
  return std::sqrt(norm_squared(value));
}
[[nodiscard]] inline bool is_finite(const Vec3 value) noexcept {
  return std::isfinite(value.x) && std::isfinite(value.y) && std::isfinite(value.z);
}
[[nodiscard]] inline Vec3 normalize(const Vec3 value) {
  if (!is_finite(value)) {
    throw std::invalid_argument("cannot normalize a non-finite vector");
  }
  const Real magnitude = norm(value);
  if (!(magnitude > 0.0) || !std::isfinite(magnitude)) {
    throw std::invalid_argument("cannot normalize a zero vector");
  }
  return value / magnitude;
}
[[nodiscard]] constexpr Vec3 project_tangent(const Vec3 value,
                                             const Vec3 unit_position) noexcept {
  return value - dot(value, unit_position) * unit_position;
}
[[nodiscard]] inline Real safe_angle(const Vec3 first, const Vec3 second) {
  const Vec3 a = normalize(first);
  const Vec3 b = normalize(second);
  return std::atan2(norm(cross(a, b)), dot(a, b));
}

}  // namespace mps
