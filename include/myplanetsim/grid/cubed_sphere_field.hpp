#pragma once

#include <cstddef>
#include <span>
#include <vector>

#include "myplanetsim/grid/cubed_sphere_grid.hpp"
#include "myplanetsim/grid/field2d.hpp"

namespace mps {

class CubedSphereField {
 public:
  CubedSphereField(Index cells_per_panel, Index halo_width = 1);

  [[nodiscard]] Index cells_per_panel() const noexcept { return n_; }
  [[nodiscard]] Index halo_width() const noexcept { return halo_; }
  [[nodiscard]] std::size_t cell_count() const noexcept;
  [[nodiscard]] Field2D<Real>& panel(Panel panel) noexcept;
  [[nodiscard]] const Field2D<Real>& panel(Panel panel) const noexcept;
  [[nodiscard]] Real& operator()(CellId cell) noexcept;
  [[nodiscard]] const Real& operator()(CellId cell) const noexcept;

  void fill(Real value);
  [[nodiscard]] std::vector<Real> flatten() const;
  void assign_flattened(std::span<const Real> values);
  void exchange_edge_halos();

 private:
  Index n_;
  Index halo_;
  std::vector<Field2D<Real>> panels_;
};

}  // namespace mps
