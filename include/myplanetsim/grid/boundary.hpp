#pragma once

#include "myplanetsim/core/types.hpp"
#include "myplanetsim/grid/field2d.hpp"

namespace mps {
namespace detail {

[[nodiscard]] constexpr Index periodic_index(const Index index,
                                             const Index extent) noexcept {
  const Index remainder = index % extent;
  return remainder < 0 ? remainder + extent : remainder;
}

}  // namespace detail

template <typename T>
void fill_periodic_halo(Field2D<T>& field) {
  const auto halo = field.halo();
  for (Index j = -halo; j < field.ny() + halo; ++j) {
    for (Index i = -halo; i < field.nx() + halo; ++i) {
      const bool is_interior = i >= 0 && i < field.nx() && j >= 0 && j < field.ny();
      if (is_interior) {
        continue;
      }

      const auto source_i = detail::periodic_index(i, field.nx());
      const auto source_j = detail::periodic_index(j, field.ny());
      field(i, j) = field(source_i, source_j);
    }
  }
}

}  // namespace mps
