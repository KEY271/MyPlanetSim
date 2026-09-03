#include <algorithm>
#include <cmath>
#include <vector>

#include "myplanetsim/numerics/compatible_operators.hpp"
#include "support/test.hpp"

namespace {

constexpr mps::Real kRadius = 6.371e6;
constexpr mps::Vec3 kRotation{0.0, 0.0, 1.0e-5};

[[nodiscard]] std::vector<mps::Vec3> solid_rotation_velocity(
    const mps::CubedSphereGrid& grid) {
  std::vector<mps::Vec3> velocity(grid.cell_count());
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
    velocity[cell] = grid.radius_m() * cross(kRotation, grid.cells()[cell].center);
  }
  return velocity;
}

struct VorticityError {
  mps::Real median;
  mps::Real maximum;
};

// Solid-rotation vorticity error. The median measures the smooth part of the
// grid; the maximum is dominated by the eight three-valent cube corners and
// their neighbours, where the dual metric stays skewed under refinement.
[[nodiscard]] VorticityError vorticity_error(const mps::Index resolution) {
  const mps::CubedSphereGrid grid(resolution, kRadius);
  const mps::CubedSphereDualTopology dual(grid);
  const auto normal_velocity =
      mps::edge_normal_velocity(grid, dual, solid_rotation_velocity(grid));
  const auto vorticity = mps::vertex_relative_vorticity(grid, dual, normal_velocity);
  std::vector<mps::Real> errors(grid.vertex_count());
  for (std::size_t vertex = 0; vertex < grid.vertex_count(); ++vertex) {
    const mps::Real exact = 2.0 * mps::dot(kRotation, grid.vertices()[vertex].position);
    errors[vertex] = std::abs(vorticity[vertex] - exact);
  }
  std::ranges::sort(errors);
  return {.median = errors[errors.size() / 2], .maximum = errors.back()};
}

[[nodiscard]] mps::Real perot_velocity_error(const mps::Index resolution) {
  const mps::CubedSphereGrid grid(resolution, kRadius);
  const mps::CubedSphereDualTopology dual(grid);
  const auto exact = solid_rotation_velocity(grid);
  const auto normal_velocity = mps::edge_normal_velocity(grid, dual, exact);
  const auto reconstructed = mps::perot_cell_velocity(grid, normal_velocity);
  mps::Real worst = 0.0;
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
    worst = std::max(worst, mps::norm(reconstructed[cell] - exact[cell]));
  }
  return worst;
}

}  // namespace

MPS_TEST_CASE("compatible edge metrics are unit vectors with bounded skewness") {
  const mps::CubedSphereGrid grid(16, kRadius);
  const mps::CubedSphereDualTopology dual(grid);
  const auto metrics = mps::compatible_edge_metrics(grid, dual);
  mps::Real smallest_normal_projection = 1.0;
  for (const auto& metric : metrics) {
    MPS_CHECK(metric.hodge_ratio > 0.0);
    MPS_CHECK_NEAR(metric.normal_projection * metric.normal_projection +
                       metric.tangent_projection * metric.tangent_projection,
                   1.0, 1.0e-12);
    smallest_normal_projection =
        std::min(smallest_normal_projection, metric.normal_projection);
  }
  // The equiangular cubed sphere is not orthogonal, so the dual edge leans away
  // from the primal edge normal by up to about twenty-five degrees.
  MPS_CHECK(smallest_normal_projection > 0.9);
  MPS_CHECK(smallest_normal_projection < 1.0);
}

MPS_TEST_CASE("Perot reconstruction recovers a solid-rotation velocity") {
  const mps::Real coarse = perot_velocity_error(8);
  const mps::Real fine = perot_velocity_error(16);
  MPS_CHECK(coarse < 0.05 * kRadius * kRotation.z);
  MPS_CHECK(fine < 0.7 * coarse);
}

MPS_TEST_CASE("dual circulation reproduces the solid-rotation vorticity") {
  const auto coarse = vorticity_error(8);
  const auto fine = vorticity_error(16);
  const mps::Real scale = 2.0 * kRotation.z;
  MPS_CHECK(coarse.median < 0.01 * scale);
  MPS_CHECK(fine.median < 0.5 * coarse.median);
  // The cube corners do not converge; they are bounded and reported instead.
  MPS_CHECK(coarse.maximum < 0.3 * scale);
  MPS_CHECK(fine.maximum < 1.2 * coarse.maximum);
}

MPS_TEST_CASE("dual circulation sums to zero over the sphere") {
  const mps::CubedSphereGrid grid(8, kRadius);
  const mps::CubedSphereDualTopology dual(grid);
  const auto normal_velocity =
      mps::edge_normal_velocity(grid, dual, solid_rotation_velocity(grid));
  const auto vorticity = mps::vertex_relative_vorticity(grid, dual, normal_velocity);
  mps::Real circulation = 0.0;
  mps::Real scale = 0.0;
  for (const auto& dual_vertex : dual.vertices()) {
    const mps::Real weighted =
        dual_vertex.area_m2 * vorticity[dual_vertex.primal_vertex];
    circulation += weighted;
    scale += std::abs(weighted);
  }
  MPS_CHECK(std::abs(circulation) < 1.0e-12 * scale);
}

MPS_TEST_CASE("resting state has planetary potential vorticity on every vertex") {
  const mps::CubedSphereGrid grid(8, kRadius);
  const mps::CubedSphereDualTopology dual(grid);
  constexpr mps::Real depth_m = 1000.0;
  const std::vector<mps::Real> depth(grid.cell_count(), depth_m);
  const std::vector<mps::Real> normal_velocity(grid.edge_count(), 0.0);
  const auto potential_vorticity =
      mps::vertex_potential_vorticity(grid, dual, depth, normal_velocity, kRotation);
  const auto vertex = mps::interpolate_cells_to_vertices(grid, dual, depth);
  for (std::size_t index = 0; index < grid.vertex_count(); ++index) {
    const mps::Real exact =
        2.0 * mps::dot(kRotation, grid.vertices()[index].position) / depth_m;
    MPS_CHECK_NEAR(potential_vorticity[index], exact, 1.0e-16);
    MPS_CHECK_NEAR(vertex[index], depth_m, 1.0e-12 * depth_m);
  }
}

MPS_TEST_CASE("edge interpolation preserves a constant potential vorticity") {
  const mps::CubedSphereGrid grid(4, kRadius);
  const std::vector<mps::Real> vertex_values(grid.vertex_count(), 3.5);
  const auto edge_values = mps::edge_potential_vorticity(grid, vertex_values);
  MPS_CHECK_EQ(edge_values.size(), grid.edge_count());
  for (const mps::Real value : edge_values) {
    MPS_CHECK_EQ(value, 3.5);
  }
}

MPS_TEST_CASE("cell kinetic energy matches the solid-rotation speed") {
  const mps::CubedSphereGrid grid(16, kRadius);
  const mps::CubedSphereDualTopology dual(grid);
  const auto velocity = solid_rotation_velocity(grid);
  const auto normal_velocity = mps::edge_normal_velocity(grid, dual, velocity);
  const auto kinetic_energy = mps::cell_kinetic_energy(grid, normal_velocity);
  const mps::Real scale = 0.5 * kRadius * kRotation.z * kRadius * kRotation.z;
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
    MPS_CHECK_NEAR(kinetic_energy[cell], 0.5 * mps::norm_squared(velocity[cell]),
                   0.02 * scale);
  }
}

MPS_TEST_CASE("compatible operators reject mismatched fields") {
  const mps::CubedSphereGrid grid(2, kRadius);
  const mps::CubedSphereDualTopology dual(grid);
  const std::vector<mps::Real> wrong(grid.edge_count() + 1, 0.0);
  MPS_CHECK_THROWS_AS(mps::vertex_relative_vorticity(grid, dual, wrong),
                      std::invalid_argument);
  MPS_CHECK_THROWS_AS(mps::perot_cell_velocity(grid, wrong), std::invalid_argument);
  MPS_CHECK_THROWS_AS(mps::edge_potential_vorticity(grid, wrong),
                      std::invalid_argument);
}

int main() { return mps::test::run_all(); }
