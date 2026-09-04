#include <cmath>
#include <vector>

#include "myplanetsim/dynamics/shallow_water_diffusion.hpp"
#include "myplanetsim/numerics/spherical_operators.hpp"
#include "support/test.hpp"

// Gates for the spherical vector Laplacian of ADR 0011,
//
//   lap(v) = grad(div v) + k_hat x grad(zeta) + v / a^2
//
// whose rotational sign and curvature term were both wrong before Phase 9. The old
// operator returned <lap(u).u>/<u.u> = +1.996 on a unit sphere for rigid rotation where
// the exact value is -1, so every check below fails against it.

namespace {

constexpr mps::Real kRadius = 3.0;
constexpr mps::Vec3 kAxis{0.0, 0.0, 1.0};

// Analytic fields built from a degree-n spherical harmonic have the eigenvalue
// (1 - n(n+1)) / a^2 for both their rotational and their divergent form.
[[nodiscard]] mps::Real eigenvalue(const mps::Real degree) {
  return (1.0 - degree * (degree + 1.0)) / (kRadius * kRadius);
}

// Surface gradient of psi_n evaluated at the unit direction r_hat, for
// psi_1 = z and psi_2 = (3 z^2 - 1) / 2 with z = k_hat . r_hat.
[[nodiscard]] mps::Vec3 harmonic_surface_gradient(const mps::Vec3 direction,
                                                  const int degree) {
  const mps::Real z = mps::dot(kAxis, direction);
  const mps::Vec3 meridional = kAxis - z * direction;
  const mps::Real factor = degree == 1 ? 1.0 : 3.0 * z;
  return (factor / kRadius) * meridional;
}

enum class Form { kRotational, kDivergent };

[[nodiscard]] mps::Vec3 harmonic_field(const mps::Vec3 direction, const int degree,
                                       const Form form) {
  const mps::Vec3 gradient = harmonic_surface_gradient(direction, degree);
  return form == Form::kDivergent ? gradient : mps::cross(direction, gradient);
}

// Area-weighted Rayleigh quotient <lap(v).v> / <v.v>. For an eigenfield this is the
// eigenvalue, and it is the quantity diffusion actually consumes.
//
// The gate is stated on this projection rather than on a pointwise norm because
// finite_volume_curl carries an O(1) error along the six panel seams (Linf ~5e-2 at
// every resolution while its interior median converges at second order). Taking a
// least-squares gradient of that seam noise makes the pointwise Laplacian error grow
// like 1/h. That is a separate defect of the curl stencil, outside the ADR 0011 scope,
// and it is recorded in the Phase 9 validation report rather than papered over here.
[[nodiscard]] mps::Real laplacian_eigenvalue(const mps::Index n, const int degree,
                                             const Form form) {
  const mps::CubedSphereGrid grid(n, kRadius);
  std::vector<mps::Vec3> field(grid.cell_count());
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
    field[cell] = harmonic_field(grid.cells()[cell].center, degree, form);
  }
  const auto laplacian = mps::finite_volume_vector_laplacian(grid, field);
  mps::Real weighted_work = 0.0;
  mps::Real weighted_energy = 0.0;
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
    const mps::Real area = grid.cells()[cell].area_m2;
    weighted_work += area * mps::dot(laplacian[cell], field[cell]);
    weighted_energy += area * mps::norm_squared(field[cell]);
  }
  return weighted_work / weighted_energy;
}

[[nodiscard]] mps::Real observed_order(const mps::Real coarse, const mps::Real fine) {
  return std::log(coarse / fine) / std::log(2.0);
}

[[nodiscard]] mps::ShallowWaterState rotational_state(const mps::CubedSphereGrid& grid,
                                                      const mps::Real depth_m) {
  mps::ShallowWaterState state;
  state.depth.assign(grid.cell_count(), depth_m);
  state.momentum.resize(grid.cell_count());
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
    state.momentum[cell] =
        depth_m * harmonic_field(grid.cells()[cell].center, 2, Form::kRotational);
  }
  return state;
}

[[nodiscard]] mps::Real kinetic_energy_rate(const mps::CubedSphereGrid& grid,
                                            const mps::DiffusionKind kind,
                                            const mps::Real coefficient) {
  constexpr mps::Real depth_m = 1000.0;
  const auto state = rotational_state(grid, depth_m);
  const auto tendency =
      mps::shallow_water_diffusion_tendency(grid, state, kind, coefficient);
  mps::Real rate = 0.0;
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
    rate += grid.cells()[cell].area_m2 *
            mps::dot(state.velocity(cell), tendency.momentum[cell]);
  }
  return rate;
}

}  // namespace

MPS_TEST_CASE("vector Laplacian reproduces the rigid-rotation eigenvalue") {
  const mps::CubedSphereGrid grid(16, kRadius);
  std::vector<mps::Vec3> velocity(grid.cell_count());
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
    velocity[cell] = kRadius * mps::cross(kAxis, grid.cells()[cell].center);
  }
  const auto laplacian = mps::finite_volume_vector_laplacian(grid, velocity);
  mps::Real work = 0.0;
  mps::Real energy = 0.0;
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
    const mps::Real area = grid.cells()[cell].area_m2;
    work += area * mps::dot(laplacian[cell], velocity[cell]);
    energy += area * mps::norm_squared(velocity[cell]);
  }
  // The pre-Phase-9 operator gave +1.996 / a^2 here.
  MPS_CHECK_NEAR(work / energy, eigenvalue(1.0), 0.02 * std::abs(eigenvalue(1.0)));
}

MPS_TEST_CASE("vector Laplacian converges to the analytic harmonic eigenvalues") {
  for (const Form form : {Form::kRotational, Form::kDivergent}) {
    for (const int degree : {1, 2}) {
      const mps::Real exact = eigenvalue(static_cast<mps::Real>(degree));
      const mps::Real coarse = std::abs(laplacian_eigenvalue(8, degree, form) - exact);
      const mps::Real medium = std::abs(laplacian_eigenvalue(16, degree, form) - exact);
      const mps::Real fine = std::abs(laplacian_eigenvalue(32, degree, form) - exact);
      MPS_CHECK(observed_order(coarse, medium) >= 1.0);
      MPS_CHECK(observed_order(medium, fine) >= 1.0);
      // The eigenvalue is negative for every degree n >= 1, so the operator can only
      // remove energy from an eigenfield. The old sign made it positive.
      MPS_CHECK(fine < 0.01 * std::abs(exact));
    }
  }
}

MPS_TEST_CASE("both diffusion kinds remove kinetic energy from a smooth field") {
  const mps::CubedSphereGrid grid(16, kRadius);
  MPS_CHECK(kinetic_energy_rate(grid, mps::DiffusionKind::kLaplacian, 1.0e-3) < 0.0);
  MPS_CHECK(kinetic_energy_rate(grid, mps::DiffusionKind::kBiharmonic, 1.0e-6) < 0.0);
}

int main() { return mps::test::run_all(); }
