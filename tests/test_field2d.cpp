#include <limits>
#include <stdexcept>

#include "myplanetsim/grid/field2d.hpp"
#include "support/test.hpp"

MPS_TEST_CASE("Field2D maps rectangular storage including halos") {
  mps::Field2D<int> field(3, 2, 1);
  MPS_CHECK_EQ(field.nx(), 3);
  MPS_CHECK_EQ(field.ny(), 2);
  MPS_CHECK_EQ(field.halo(), 1);
  MPS_CHECK_EQ(field.storage_nx(), 5);
  MPS_CHECK_EQ(field.storage_ny(), 4);
  MPS_CHECK_EQ(field.size(), 20U);

  for (mps::Index j = -1; j <= 2; ++j) {
    for (mps::Index i = -1; i <= 3; ++i) {
      field.at(i, j) = static_cast<int>(100 * j + i);
    }
  }
  for (mps::Index j = -1; j <= 2; ++j) {
    for (mps::Index i = -1; i <= 3; ++i) {
      MPS_CHECK_EQ(field.at(i, j), static_cast<int>(100 * j + i));
    }
  }
}

MPS_TEST_CASE("Field2D exposes storage and interior rows") {
  mps::Field2D<double> field(4, 3, 2);
  field.fill(7.0);
  MPS_CHECK_EQ(field.storage().size(), field.size());
  for (mps::Index j = 0; j < field.ny(); ++j) {
    MPS_CHECK_EQ(field.interior_row(j).size(), 4U);
    for (const auto value : field.interior_row(j)) {
      MPS_CHECK_EQ(value, 7.0);
    }
  }
}

MPS_TEST_CASE("Field2D copies deeply and moves ownership") {
  mps::Field2D<int> original(2, 2);
  original.fill(3);
  auto copy = original;
  copy.at(0, 0) = 9;
  MPS_CHECK_EQ(original.at(0, 0), 3);
  MPS_CHECK_EQ(copy.at(0, 0), 9);

  auto moved = std::move(copy);
  MPS_CHECK_EQ(moved.at(0, 0), 9);
}

MPS_TEST_CASE("Field2D rejects invalid extents and indices") {
  MPS_CHECK_THROWS_AS((mps::Field2D<int>(0, 1)), std::invalid_argument);
  MPS_CHECK_THROWS_AS((mps::Field2D<int>(1, -1)), std::invalid_argument);
  MPS_CHECK_THROWS_AS((mps::Field2D<int>(1, 1, -1)), std::invalid_argument);

  mps::Field2D<int> field(2, 3, 1);
  MPS_CHECK_THROWS_AS(field.at(-2, 0), std::out_of_range);
  MPS_CHECK_THROWS_AS(field.at(3, 0), std::out_of_range);
  MPS_CHECK_THROWS_AS(field.at(0, -2), std::out_of_range);
  MPS_CHECK_THROWS_AS(field.at(0, 4), std::out_of_range);
  MPS_CHECK_THROWS_AS(field.interior_row(3), std::out_of_range);

  const auto maximum = std::numeric_limits<mps::Index>::max();
  MPS_CHECK_THROWS_AS((mps::Field2D<int>(maximum, 1, 1)), std::length_error);
}

int main() { return mps::test::run_all(); }
