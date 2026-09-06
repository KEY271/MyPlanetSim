#pragma once

#include <span>
#include <vector>

#include "myplanetsim/grid/cubed_sphere_grid.hpp"
#include "myplanetsim/numerics/gmres.hpp"

namespace mps {

// Cell-centred finite-volume H = I - coefficient * Laplacian. The edge gradient and
// cell divergence share one oriented flux, so constants are preserved and the
// area-weighted global Laplacian sum is zero to roundoff.
class FiniteVolumeHelmholtzOperator {
 public:
  FiniteVolumeHelmholtzOperator(const CubedSphereGrid& grid, Real coefficient_m2);

  void apply(std::span<const Real> input, std::span<Real> output) const;
  void apply_laplacian(std::span<const Real> input, std::span<Real> output) const;
  void apply_jacobi_preconditioner(std::span<const Real> input,
                                   std::span<Real> output) const;

  [[nodiscard]] Real coefficient_m2() const noexcept { return coefficient_m2_; }
  [[nodiscard]] std::span<const Real> inverse_diagonal() const noexcept {
    return inverse_diagonal_;
  }

 private:
  const CubedSphereGrid* grid_;
  Real coefficient_m2_;
  std::vector<Real> inverse_diagonal_;
};

[[nodiscard]] GmresResult solve_finite_volume_helmholtz(
    const FiniteVolumeHelmholtzOperator& helmholtz,
    std::span<const Real> right_hand_side, std::span<Real> solution,
    const GmresOptions& options, GmresWorkspace& workspace);

}  // namespace mps
