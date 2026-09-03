#include <cmath>
#include <vector>

#include "myplanetsim/grid/cubed_sphere_field.hpp"
#include "support/test.hpp"

MPS_TEST_CASE("panel field flattening is panel major and reversible") {
  mps::CubedSphereField field(4, 2);
  std::vector<double> expected(field.cell_count());
  for (std::size_t index = 0; index < expected.size(); ++index) {
    expected[index] = static_cast<double>(index) + 0.25;
  }
  field.assign_flattened(expected);
  MPS_CHECK(field.flatten() == expected);
  auto copy = field;
  copy({mps::Panel::kPositiveX, 0, 0}) = -1.0;
  MPS_CHECK(field({mps::Panel::kPositiveX, 0, 0}) !=
            copy({mps::Panel::kPositiveX, 0, 0}));
  MPS_CHECK_THROWS_AS(field.assign_flattened(std::vector<double>(3)),
                      std::invalid_argument);
}

MPS_TEST_CASE("edge halos are exchanged and corners stay poisoned") {
  mps::CubedSphereField field(5, 2);
  for (const auto panel : mps::kPanels) {
    for (mps::Index j = 0; j < 5; ++j) {
      for (mps::Index i = 0; i < 5; ++i) {
        field({panel, i, j}) = 100.0 * static_cast<double>(mps::panel_index(panel)) +
                               10.0 * static_cast<double>(j) + static_cast<double>(i);
      }
    }
  }
  const auto interior = field.flatten();
  field.exchange_edge_halos();
  MPS_CHECK(field.flatten() == interior);
  for (const auto panel : mps::kPanels) {
    const auto& values = field.panel(panel);
    for (mps::Index layer = 1; layer <= 2; ++layer) {
      for (mps::Index along = 0; along < 5; ++along) {
        MPS_CHECK(std::isfinite(values(-layer, along)));
        MPS_CHECK(std::isfinite(values(4 + layer, along)));
        MPS_CHECK(std::isfinite(values(along, -layer)));
        MPS_CHECK(std::isfinite(values(along, 4 + layer)));
      }
    }
    MPS_CHECK(std::isnan(values(-1, -1)));
    MPS_CHECK(std::isnan(values(5, 5)));
  }
}

int main() { return mps::test::run_all(); }
