#include "myplanetsim/numerics/spherical_operators.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "myplanetsim/geometry/cubed_sphere.hpp"

namespace mps {
namespace {

void validate_cell_values(const CubedSphereGrid& grid,
                          const std::span<const Real> values) {
  if (values.size() != grid.cell_count()) {
    throw std::invalid_argument("cell field size does not match grid");
  }
  if (!std::ranges::all_of(values,
                           [](const Real value) { return std::isfinite(value); })) {
    throw std::invalid_argument("cell field contains a non-finite value");
  }
}

[[nodiscard]] Vec3 logarithmic_displacement(const Vec3 from, const Vec3 to,
                                            const Real radius_m) {
  const Real cosine = std::clamp(dot(from, to), -1.0, 1.0);
  const Real angle = std::acos(cosine);
  const Real sine = std::sin(angle);
  if (!(sine > 0.0)) {
    throw std::runtime_error("degenerate neighbor displacement");
  }
  return (radius_m * angle / sine) * (to - cosine * from);
}

[[nodiscard]] Vec3 parallel_transport(const Vec3 vector, const Vec3 from,
                                      const Vec3 to) {
  const Real denominator = 1.0 + dot(from, to);
  if (!(denominator > 0.0)) {
    throw std::runtime_error("cannot parallel transport between antipodal points");
  }
  return vector - (dot(vector, to) / denominator) * (from + to);
}

void validate_cell_vectors(const CubedSphereGrid& grid,
                           const std::span<const Vec3> values) {
  if (values.size() != grid.cell_count()) {
    throw std::invalid_argument("vector field size does not match grid");
  }
  if (!std::ranges::all_of(values, is_finite)) {
    throw std::invalid_argument("vector field contains a non-finite value");
  }
}

}  // namespace

void scatter_oriented_edge_flux(const CubedSphereGrid& grid,
                                const std::span<const Real> edge_flux,
                                const std::span<Real> cell_tendency) {
  if (edge_flux.size() != grid.edge_count() ||
      cell_tendency.size() != grid.cell_count()) {
    throw std::invalid_argument("flux accumulator size does not match grid");
  }
  std::fill(cell_tendency.begin(), cell_tendency.end(), 0.0);
  for (const auto& edge : grid.edges()) {
    const Real flux = edge_flux[edge.id];
    if (!std::isfinite(flux)) {
      throw std::invalid_argument("edge flux contains a non-finite value");
    }
    cell_tendency[grid.cell_index(edge.left_cell)] -= flux;
    cell_tendency[grid.cell_index(edge.right_cell)] += flux;
  }
}

std::vector<Real> finite_volume_divergence(
    const CubedSphereGrid& grid, const std::span<const Real> oriented_edge_flux) {
  std::vector<Real> negative_outflow(grid.cell_count());
  scatter_oriented_edge_flux(grid, oriented_edge_flux, negative_outflow);
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
    negative_outflow[cell] = -negative_outflow[cell] / grid.cells()[cell].area_m2;
  }
  return negative_outflow;
}

std::vector<Vec3> least_squares_gradient(const CubedSphereGrid& grid,
                                         const std::span<const Real> cell_values) {
  validate_cell_values(grid, cell_values);
  std::vector<Vec3> gradients(grid.cell_count());
  for (std::size_t index = 0; index < grid.cell_count(); ++index) {
    const auto& cell = grid.cells()[index];
    const auto coordinates = inverse_map(cell.center);
    const auto basis =
        tangent_basis(coordinates.panel, coordinates.alpha, coordinates.beta);
    Real a00 = 0.0;
    Real a01 = 0.0;
    Real a11 = 0.0;
    Real b0 = 0.0;
    Real b1 = 0.0;
    for (const std::size_t edge_id : grid.cell_edges(cell.id)) {
      const CellId neighbor = grid.neighbor_across(edge_id, cell.id);
      const Vec3 displacement = logarithmic_displacement(
          cell.center, grid.cell(neighbor).center, grid.radius_m());
      const Real x = dot(displacement, basis.alpha);
      const Real y = dot(displacement, basis.beta);
      const Real difference =
          cell_values[grid.cell_index(neighbor)] - cell_values[index];
      const Real weight = 1.0 / std::max(x * x + y * y, 1.0e-300);
      a00 += weight * x * x;
      a01 += weight * x * y;
      a11 += weight * y * y;
      b0 += weight * x * difference;
      b1 += weight * y * difference;
    }
    const Real determinant = a00 * a11 - a01 * a01;
    if (!(determinant > 1.0e-12 * a00 * a11)) {
      throw std::runtime_error("least-squares gradient stencil is rank deficient");
    }
    const Real gx = (a11 * b0 - a01 * b1) / determinant;
    const Real gy = (a00 * b1 - a01 * b0) / determinant;
    gradients[index] = gx * basis.alpha + gy * basis.beta;
  }
  return gradients;
}

std::vector<Real> finite_volume_laplacian(const CubedSphereGrid& grid,
                                          const std::span<const Real> cell_values) {
  validate_cell_values(grid, cell_values);
  std::vector<Real> edge_flux(grid.edge_count());
  for (const auto& edge : grid.edges()) {
    const auto left = grid.cell_index(edge.left_cell);
    const auto right = grid.cell_index(edge.right_cell);
    const Real distance = grid.radius_m() * safe_angle(grid.cells()[left].center,
                                                       grid.cells()[right].center);
    edge_flux[edge.id] =
        (cell_values[right] - cell_values[left]) * edge.length_m / distance;
  }
  return finite_volume_divergence(grid, edge_flux);
}

std::vector<Real> finite_volume_curl(const CubedSphereGrid& grid,
                                     const std::span<const Vec3> cell_vectors) {
  validate_cell_vectors(grid, cell_vectors);
  std::vector<Real> circulation(grid.cell_count(), 0.0);
  for (const auto& edge : grid.edges()) {
    const auto left = grid.cell_index(edge.left_cell);
    const auto right = grid.cell_index(edge.right_cell);
    const Vec3 tangent =
        normalize(project_tangent(grid.vertices()[edge.second_vertex].position -
                                      grid.vertices()[edge.first_vertex].position,
                                  edge.center));
    const Vec3 average = 0.5 * (cell_vectors[left] + cell_vectors[right]);
    const Real integral = dot(average, tangent) * edge.length_m;
    circulation[left] += integral;
    circulation[right] -= integral;
  }
  for (std::size_t cell = 0; cell < circulation.size(); ++cell) {
    circulation[cell] /= grid.cells()[cell].area_m2;
  }
  return circulation;
}

std::vector<Real> finite_volume_vector_divergence(
    const CubedSphereGrid& grid, const std::span<const Vec3> cell_vectors) {
  validate_cell_vectors(grid, cell_vectors);
  std::vector<Real> edge_flux(grid.edge_count());
  for (const auto& edge : grid.edges()) {
    const std::size_t left = grid.cell_index(edge.left_cell);
    const std::size_t right = grid.cell_index(edge.right_cell);
    const Vec3 average = 0.5 * (cell_vectors[left] + cell_vectors[right]);
    edge_flux[edge.id] = dot(average, edge.outward_normal_from_left) * edge.length_m;
  }
  return finite_volume_divergence(grid, edge_flux);
}

std::vector<TangentVectorGradient> least_squares_vector_gradient(
    const CubedSphereGrid& grid, const std::span<const Vec3> cell_vectors) {
  validate_cell_vectors(grid, cell_vectors);
  std::vector<TangentVectorGradient> gradients(grid.cell_count());
  for (std::size_t index = 0; index < grid.cell_count(); ++index) {
    const auto& cell = grid.cells()[index];
    const auto coordinates = inverse_map(cell.center);
    const auto basis =
        tangent_basis(coordinates.panel, coordinates.alpha, coordinates.beta);
    Real a00 = 0.0;
    Real a01 = 0.0;
    Real a11 = 0.0;
    Vec3 b0{};
    Vec3 b1{};
    for (const std::size_t edge_id : grid.cell_edges(cell.id)) {
      const CellId neighbor_id = grid.neighbor_across(edge_id, cell.id);
      const std::size_t neighbor = grid.cell_index(neighbor_id);
      const Vec3 displacement = logarithmic_displacement(
          cell.center, grid.cells()[neighbor].center, grid.radius_m());
      const Real x = dot(displacement, basis.alpha);
      const Real y = dot(displacement, basis.beta);
      const Vec3 transported = parallel_transport(
          project_tangent(cell_vectors[neighbor], grid.cells()[neighbor].center),
          grid.cells()[neighbor].center, cell.center);
      const Vec3 difference =
          transported - project_tangent(cell_vectors[index], cell.center);
      const Real weight = 1.0 / std::max(x * x + y * y, 1.0e-300);
      a00 += weight * x * x;
      a01 += weight * x * y;
      a11 += weight * y * y;
      b0 = b0 + weight * x * difference;
      b1 = b1 + weight * y * difference;
    }
    const Real determinant = a00 * a11 - a01 * a01;
    if (!(determinant > 1.0e-12 * a00 * a11)) {
      throw std::runtime_error("vector gradient stencil is rank deficient");
    }
    gradients[index] = {
        .alpha_derivative = (a11 * b0 - a01 * b1) / determinant,
        .beta_derivative = (a00 * b1 - a01 * b0) / determinant,
    };
  }
  return gradients;
}

Vec3 reconstruct_tangent_vector(const CubedSphereGrid& grid,
                                const std::size_t cell_index, const Vec3 cell_value,
                                const Vec3 face_position,
                                const TangentVectorGradient& gradient) {
  if (cell_index >= grid.cell_count()) {
    throw std::out_of_range("vector reconstruction cell is out of range");
  }
  const auto& cell = grid.cells()[cell_index];
  const auto coordinates = inverse_map(cell.center);
  const auto basis =
      tangent_basis(coordinates.panel, coordinates.alpha, coordinates.beta);
  const Vec3 face = normalize(face_position);
  const Vec3 displacement =
      logarithmic_displacement(cell.center, face, grid.radius_m());
  if (!is_finite(cell_value)) {
    throw std::invalid_argument("vector reconstruction value is non-finite");
  }
  const Vec3 increment = project_tangent(gradient.alpha_derivative, cell.center) *
                             dot(displacement, basis.alpha) +
                         project_tangent(gradient.beta_derivative, cell.center) *
                             dot(displacement, basis.beta);
  const Vec3 reconstructed = parallel_transport(
      project_tangent(cell_value, cell.center) + increment, cell.center, face);
  return project_tangent(reconstructed, face);
}

EdgeTangentBasis edge_tangent_basis(const EdgeGeometry& edge) {
  const Vec3 normal =
      normalize(project_tangent(edge.outward_normal_from_left, edge.center));
  return {.normal = normal, .tangent = normalize(cross(edge.center, normal))};
}

std::vector<Real> shallow_water_potential_vorticity(
    const CubedSphereGrid& grid, const std::span<const Real> depth,
    const std::span<const Vec3> velocity, const Vec3 rotation_vector_rad_s) {
  validate_cell_values(grid, depth);
  validate_cell_vectors(grid, velocity);
  if (!is_finite(rotation_vector_rad_s)) {
    throw std::invalid_argument("rotation vector is non-finite");
  }
  auto pv = finite_volume_curl(grid, velocity);
  for (std::size_t cell = 0; cell < pv.size(); ++cell) {
    if (!(depth[cell] > 0.0)) {
      throw std::invalid_argument("potential vorticity requires positive depth");
    }
    pv[cell] =
        (pv[cell] + 2.0 * dot(rotation_vector_rad_s, grid.cells()[cell].center)) /
        depth[cell];
  }
  return pv;
}

std::vector<Vec3> finite_volume_vector_laplacian(
    const CubedSphereGrid& grid, const std::span<const Vec3> cell_vectors) {
  const auto divergence = finite_volume_vector_divergence(grid, cell_vectors);
  const auto vorticity = finite_volume_curl(grid, cell_vectors);
  const auto divergence_gradient = least_squares_gradient(grid, divergence);
  const auto vorticity_gradient = least_squares_gradient(grid, vorticity);
  std::vector<Vec3> result(grid.cell_count());
  for (std::size_t cell = 0; cell < result.size(); ++cell) {
    result[cell] =
        project_tangent(divergence_gradient[cell] -
                            cross(grid.cells()[cell].center, vorticity_gradient[cell]),
                        grid.cells()[cell].center);
  }
  return result;
}

}  // namespace mps
