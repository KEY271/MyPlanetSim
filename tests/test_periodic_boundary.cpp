#include "myplanetsim/grid/boundary.hpp"
#include "myplanetsim/grid/field2d.hpp"
#include "support/test.hpp"

namespace {

[[nodiscard]] int interior_value(const mps::Index i, const mps::Index j) {
  return static_cast<int>(100 * j + i);
}

void check_periodic_field(const mps::Index nx, const mps::Index ny,
                          const mps::Index halo) {
  mps::Field2D<int> field(nx, ny, halo);
  field.fill(-1);
  for (mps::Index j = 0; j < ny; ++j) {
    for (mps::Index i = 0; i < nx; ++i) {
      field(i, j) = interior_value(i, j);
    }
  }

  mps::fill_periodic_halo(field);

  for (mps::Index j = -halo; j < ny + halo; ++j) {
    for (mps::Index i = -halo; i < nx + halo; ++i) {
      const auto source_i = mps::detail::periodic_index(i, nx);
      const auto source_j = mps::detail::periodic_index(j, ny);
      MPS_CHECK_EQ(field(i, j), interior_value(source_i, source_j));
    }
  }
}

}  // namespace

MPS_TEST_CASE("periodic halo fills rectangular edges and corners") {
  check_periodic_field(4, 3, 1);
  check_periodic_field(5, 2, 2);
}

MPS_TEST_CASE("periodic halo supports widths larger than an extent") {
  check_periodic_field(2, 3, 4);
}

MPS_TEST_CASE("periodic halo preserves constant fields") {
  mps::Field2D<double> field(3, 7, 2);
  field.fill(4.5);
  mps::fill_periodic_halo(field);
  for (const auto value : field.storage()) {
    MPS_CHECK_EQ(value, 4.5);
  }
}

int main() { return mps::test::run_all(); }
