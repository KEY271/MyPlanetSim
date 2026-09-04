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
    const auto& cached = grid.edge_cache()[edge.id];
    const Real flux = edge_flux[edge.id];
    if (!std::isfinite(flux)) {
      throw std::invalid_argument("edge flux contains a non-finite value");
    }
    cell_tendency[cached.left_cell] -= flux;
    cell_tendency[cached.right_cell] += flux;
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
    const auto& cached_cell = grid.cell_cache()[index];
    TangentComponents components{};
    for (const auto& cached_edge : cached_cell.edges) {
      const Real difference = cell_values[cached_edge.neighbor] - cell_values[index];
      components.alpha += cached_edge.least_squares_weight_m_inverse.alpha * difference;
      components.beta += cached_edge.least_squares_weight_m_inverse.beta * difference;
    }
    gradients[index] = components.alpha * cached_cell.basis.alpha +
                       components.beta * cached_cell.basis.beta;
  }
  return gradients;
}

std::vector<Real> finite_volume_laplacian(const CubedSphereGrid& grid,
                                          const std::span<const Real> cell_values) {
  validate_cell_values(grid, cell_values);
  std::vector<Real> edge_flux(grid.edge_count());
  for (const auto& edge : grid.edges()) {
    const auto& cached = grid.edge_cache()[edge.id];
    const auto left = cached.left_cell;
    const auto right = cached.right_cell;
    edge_flux[edge.id] = (cell_values[right] - cell_values[left]) * edge.length_m /
                         cached.center_distance_m;
  }
  return finite_volume_divergence(grid, edge_flux);
}

std::vector<Real> finite_volume_curl(const CubedSphereGrid& grid,
                                     const std::span<const Vec3> cell_vectors) {
  validate_cell_vectors(grid, cell_vectors);
  std::vector<Real> circulation(grid.cell_count(), 0.0);
  for (const auto& edge : grid.edges()) {
    const auto& cached = grid.edge_cache()[edge.id];
    const auto left = cached.left_cell;
    const auto right = cached.right_cell;
    const Vec3 average = 0.5 * (cell_vectors[left] + cell_vectors[right]);
    const Real integral = dot(average, cached.circulation_tangent) * edge.length_m;
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
    const auto& cached = grid.edge_cache()[edge.id];
    const std::size_t left = cached.left_cell;
    const std::size_t right = cached.right_cell;
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
    const auto& cached_cell = grid.cell_cache()[index];
    TangentVectorGradient gradient{};
    for (const auto& cached_edge : cached_cell.edges) {
      const std::size_t neighbor = cached_edge.neighbor;
      const Vec3 transported = parallel_transport(
          project_tangent(cell_vectors[neighbor], grid.cells()[neighbor].center),
          grid.cells()[neighbor].center, cell.center);
      const Vec3 difference =
          transported - project_tangent(cell_vectors[index], cell.center);
      gradient.alpha_derivative =
          gradient.alpha_derivative +
          cached_edge.least_squares_weight_m_inverse.alpha * difference;
      gradient.beta_derivative =
          gradient.beta_derivative +
          cached_edge.least_squares_weight_m_inverse.beta * difference;
    }
    gradients[index] = gradient;
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
  const Vec3 face = normalize(face_position);
  const Vec3 displacement =
      logarithmic_displacement(cell.center, face, grid.radius_m());
  return reconstruct_tangent_vector_cached(grid, cell_index, cell_value, face,
                                           displacement, gradient);
}

Vec3 reconstruct_tangent_vector_cached(const CubedSphereGrid& grid,
                                       const std::size_t cell_index,
                                       const Vec3 cell_value, const Vec3 face_position,
                                       const Vec3 face_displacement_m,
                                       const TangentVectorGradient& gradient) {
  if (cell_index >= grid.cell_count()) {
    throw std::out_of_range("vector reconstruction cell is out of range");
  }
  const auto& cell = grid.cells()[cell_index];
  const auto& basis = grid.cell_cache()[cell_index].basis;
  const Vec3 face = normalize(face_position);
  if (!is_finite(cell_value)) {
    throw std::invalid_argument("vector reconstruction value is non-finite");
  }
  const Vec3 increment = project_tangent(gradient.alpha_derivative, cell.center) *
                             dot(face_displacement_m, basis.alpha) +
                         project_tangent(gradient.beta_derivative, cell.center) *
                             dot(face_displacement_m, basis.beta);
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
  // lap(v) = grad(div v) + k_hat x grad(zeta) + v / a^2 (ADR 0011). The curvature term
  // is the same +v/a^2 for the rotational and the divergent part.
  const Real inverse_radius_squared = 1.0 / (grid.radius_m() * grid.radius_m());
  std::vector<Vec3> result(grid.cell_count());
  for (std::size_t cell = 0; cell < result.size(); ++cell) {
    const Vec3 center = grid.cells()[cell].center;
    result[cell] = project_tangent(
        divergence_gradient[cell] + cross(center, vorticity_gradient[cell]) +
            inverse_radius_squared * project_tangent(cell_vectors[cell], center),
        center);
  }
  return result;
}

}  // namespace mps
