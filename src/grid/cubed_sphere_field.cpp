#include "myplanetsim/grid/cubed_sphere_field.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>

namespace mps {

CubedSphereField::CubedSphereField(const Index cells_per_panel, const Index halo_width)
    : n_(cells_per_panel), halo_(halo_width) {
  if (n_ <= 0) {
    throw std::invalid_argument("cells_per_panel must be positive");
  }
  if (halo_ < 0 || halo_ > n_) {
    throw std::invalid_argument("halo width must be in [0, cells_per_panel]");
  }
  panels_.reserve(kPanels.size());
  for ([[maybe_unused]] const Panel panel_value : kPanels) {
    panels_.emplace_back(n_, n_, halo_);
  }
  fill(0.0);
}

std::size_t CubedSphereField::cell_count() const noexcept {
  const auto n = static_cast<std::size_t>(n_);
  return kPanels.size() * n * n;
}

Field2D<Real>& CubedSphereField::panel(const Panel panel_value) noexcept {
  return panels_[panel_index(panel_value)];
}

const Field2D<Real>& CubedSphereField::panel(const Panel panel_value) const noexcept {
  return panels_[panel_index(panel_value)];
}

Real& CubedSphereField::operator()(const CellId cell) noexcept {
  return panel(cell.panel)(cell.i, cell.j);
}

const Real& CubedSphereField::operator()(const CellId cell) const noexcept {
  return panel(cell.panel)(cell.i, cell.j);
}

void CubedSphereField::fill(const Real value) {
  for (auto& panel_field : panels_) {
    panel_field.fill(value);
  }
}

std::vector<Real> CubedSphereField::flatten() const {
  std::vector<Real> values;
  values.reserve(cell_count());
  for (const Panel panel_value : kPanels) {
    const auto& panel_field = panel(panel_value);
    for (Index j = 0; j < n_; ++j) {
      const auto row = panel_field.interior_row(j);
      values.insert(values.end(), row.begin(), row.end());
    }
  }
  return values;
}

void CubedSphereField::assign_flattened(const std::span<const Real> values) {
  if (values.size() != cell_count()) {
    throw std::invalid_argument("flattened cubed-sphere field has wrong size");
  }
  std::size_t offset = 0;
  for (const Panel panel_value : kPanels) {
    auto& panel_field = panel(panel_value);
    for (Index j = 0; j < n_; ++j) {
      auto row = panel_field.interior_row(j);
      std::copy_n(values.begin() + static_cast<std::ptrdiff_t>(offset), row.size(),
                  row.begin());
      offset += row.size();
    }
  }
}

void CubedSphereField::exchange_edge_halos() {
  if (halo_ == 0) {
    return;
  }
  const Real poison = std::numeric_limits<Real>::quiet_NaN();
  for (auto& panel_field : panels_) {
    for (Index j = -halo_; j < n_ + halo_; ++j) {
      for (Index i = -halo_; i < n_ + halo_; ++i) {
        if ((i < 0 || i >= n_) && (j < 0 || j >= n_)) {
          panel_field(i, j) = poison;
        }
      }
    }
  }

  const Real delta = std::numbers::pi_v<Real> / (2.0 * static_cast<Real>(n_));
  for (const Panel source_panel : kPanels) {
    auto& destination = panel(source_panel);
    for (Index layer = 1; layer <= halo_; ++layer) {
      for (Index along = 0; along < n_; ++along) {
        const std::array<std::pair<Index, Index>, 4> ghosts{{{-layer, along},
                                                             {n_ - 1 + layer, along},
                                                             {along, -layer},
                                                             {along, n_ - 1 + layer}}};
        for (const auto [i, j] : ghosts) {
          const Real alpha =
              -std::numbers::pi_v<Real> / 4.0 + (static_cast<Real>(i) + 0.5) * delta;
          const Real beta =
              -std::numbers::pi_v<Real> / 4.0 + (static_cast<Real>(j) + 0.5) * delta;
          const auto mapped =
              inverse_map(map_to_unit_sphere(source_panel, alpha, beta));
          Index source_i = static_cast<Index>(
              std::floor((mapped.alpha + std::numbers::pi_v<Real> / 4.0) / delta));
          Index source_j = static_cast<Index>(
              std::floor((mapped.beta + std::numbers::pi_v<Real> / 4.0) / delta));
          source_i = std::clamp(source_i, Index{0}, n_ - 1);
          source_j = std::clamp(source_j, Index{0}, n_ - 1);
          destination(i, j) = panel(mapped.panel)(source_i, source_j);
        }
      }
    }
  }
}

}  // namespace mps
