#include <vector>

#include "myplanetsim/core/field_views.hpp"
#include "support/test.hpp"

MPS_TEST_CASE("cell-column field and tracer views preserve canonical order") {
  std::vector<double> values(2 * 3 * 5);
  for (std::size_t i = 0; i < values.size(); ++i) values[i] = static_cast<double>(i);
  auto field = mps::make_cell_column_field_view<double>(values, 2, 3, 5);

  MPS_CHECK_EQ(field(1, 2, 4), values[(1 * 3 + 2) * 5 + 4]);
  auto column = field.column(1, 2);
  MPS_CHECK_EQ(column.size(), 5U);
  MPS_CHECK_EQ(column.stride(), 1);
  column[3] = -7.0;
  MPS_CHECK_EQ(values[(1 * 3 + 2) * 5 + 3], -7.0);
  MPS_CHECK_THROWS_AS(field.at(2, 0, 0), std::out_of_range);
}

MPS_TEST_CASE("blocked field maps tail lanes without changing logical values") {
  constexpr std::size_t components = 2;
  constexpr std::size_t cells = 5;
  constexpr std::size_t levels = 3;
  for (const std::size_t lanes : {4U, 8U, 16U}) {
    const auto blocks = (cells + lanes - 1) / lanes;
    std::vector<double> storage(components * blocks * levels * lanes, -1.0);
    mps::BlockedField3DView<double> blocked(storage, components, cells, levels, lanes);
    for (std::size_t component = 0; component < components; ++component)
      for (std::size_t cell = 0; cell < cells; ++cell)
        for (std::size_t level = 0; level < levels; ++level)
          blocked(component, cell, level) =
              static_cast<double>((component * cells + cell) * levels + level);
    for (std::size_t component = 0; component < components; ++component)
      for (std::size_t cell = 0; cell < cells; ++cell)
        for (std::size_t level = 0; level < levels; ++level)
          MPS_CHECK_EQ(
              blocked(component, cell, level),
              static_cast<double>((component * cells + cell) * levels + level));
    MPS_CHECK_EQ(blocked.padded_cells(), blocks * lanes);
  }
}

int main() { return mps::test::run_all(); }
