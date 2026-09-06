#include "myplanetsim/dynamics/dry_hydrostatic_fast_modes.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <numeric>
#include <ostream>
#include <stdexcept>

#include "myplanetsim/core/validation.hpp"
#include "myplanetsim/dynamics/dry_hydrostatic_fast_operator.hpp"
#include "myplanetsim/vertical/vertical_transport.hpp"

namespace mps {
namespace {

[[nodiscard]] std::size_t matrix_offset(const std::size_t row, const std::size_t column,
                                        const std::size_t size) {
  return row * size + column;
}

[[nodiscard]] std::vector<Real> identity_matrix(const std::size_t size) {
  std::vector<Real> result(size * size, 0.0);
  for (std::size_t index = 0; index < size; ++index)
    result[matrix_offset(index, index, size)] = 1.0;
  return result;
}

[[nodiscard]] std::vector<Real> multiply_matrices(const std::span<const Real> left,
                                                  const std::span<const Real> right,
                                                  const std::size_t size) {
  std::vector<Real> result(size * size, 0.0);
  for (std::size_t row = 0; row < size; ++row) {
    for (std::size_t inner = 0; inner < size; ++inner) {
      const Real value = left[matrix_offset(row, inner, size)];
      for (std::size_t column = 0; column < size; ++column)
        result[matrix_offset(row, column, size)] +=
            value * right[matrix_offset(inner, column, size)];
    }
  }
  return result;
}

// Givens QR is sufficient here because the vertical matrices are at most a few
// dozen levels and are constructed only once by the driver.
void qr_factorize(const std::span<const Real> matrix, const std::size_t size,
                  std::vector<Real>& q, std::vector<Real>& r) {
  r.assign(matrix.begin(), matrix.end());
  q = identity_matrix(size);
  for (std::size_t column = 0; column < size; ++column) {
    for (std::size_t reverse = 0; reverse + column + 1 < size; ++reverse) {
      const std::size_t lower = size - reverse - 1;
      if (lower <= column) break;
      const std::size_t upper = lower - 1;
      const Real a = r[matrix_offset(upper, column, size)];
      const Real b = r[matrix_offset(lower, column, size)];
      const Real magnitude = std::hypot(a, b);
      if (magnitude == 0.0) continue;
      const Real cosine = a / magnitude;
      const Real sine = b / magnitude;
      for (std::size_t j = column; j < size; ++j) {
        const Real first = r[matrix_offset(upper, j, size)];
        const Real second = r[matrix_offset(lower, j, size)];
        r[matrix_offset(upper, j, size)] = cosine * first + sine * second;
        r[matrix_offset(lower, j, size)] = -sine * first + cosine * second;
      }
      for (std::size_t i = 0; i < size; ++i) {
        const Real first = q[matrix_offset(i, upper, size)];
        const Real second = q[matrix_offset(i, lower, size)];
        q[matrix_offset(i, upper, size)] = cosine * first + sine * second;
        q[matrix_offset(i, lower, size)] = -sine * first + cosine * second;
      }
    }
  }
}

[[nodiscard]] std::vector<Real> invert_matrix(std::vector<Real> matrix,
                                              const std::size_t size) {
  auto inverse = identity_matrix(size);
  for (std::size_t column = 0; column < size; ++column) {
    std::size_t pivot = column;
    for (std::size_t row = column + 1; row < size; ++row) {
      if (std::abs(matrix[matrix_offset(row, column, size)]) >
          std::abs(matrix[matrix_offset(pivot, column, size)]))
        pivot = row;
    }
    const Real diagonal = matrix[matrix_offset(pivot, column, size)];
    if (!(std::abs(diagonal) > 64.0 * std::numeric_limits<Real>::epsilon()) ||
        !std::isfinite(diagonal))
      throw std::runtime_error("dry vertical eigenvector matrix is singular");
    if (pivot != column) {
      for (std::size_t j = 0; j < size; ++j) {
        std::swap(matrix[matrix_offset(column, j, size)],
                  matrix[matrix_offset(pivot, j, size)]);
        std::swap(inverse[matrix_offset(column, j, size)],
                  inverse[matrix_offset(pivot, j, size)]);
      }
    }
    const Real scale = matrix[matrix_offset(column, column, size)];
    for (std::size_t j = 0; j < size; ++j) {
      matrix[matrix_offset(column, j, size)] /= scale;
      inverse[matrix_offset(column, j, size)] /= scale;
    }
    for (std::size_t row = 0; row < size; ++row) {
      if (row == column) continue;
      const Real factor = matrix[matrix_offset(row, column, size)];
      for (std::size_t j = 0; j < size; ++j) {
        matrix[matrix_offset(row, j, size)] -=
            factor * matrix[matrix_offset(column, j, size)];
        inverse[matrix_offset(row, j, size)] -=
            factor * inverse[matrix_offset(column, j, size)];
      }
    }
  }
  return inverse;
}

struct RealEigenDecomposition {
  std::vector<Real> eigenvalues;
  std::vector<Real> right_eigenvectors;
  std::vector<Real> inverse_eigenvectors;
};

void reduce_to_upper_hessenberg(std::vector<Real>& matrix, std::vector<Real>& transform,
                                const std::size_t size) {
  std::vector<Real> reflector(size, 0.0);
  for (std::size_t column = 0; column + 2 < size; ++column) {
    Real vector_norm = 0.0;
    for (std::size_t row = column + 1; row < size; ++row)
      vector_norm = std::hypot(vector_norm, matrix[matrix_offset(row, column, size)]);
    if (vector_norm == 0.0) continue;
    std::fill(reflector.begin(), reflector.end(), 0.0);
    const Real first = matrix[matrix_offset(column + 1, column, size)];
    reflector[column + 1] =
        first + std::copysign(vector_norm, first == 0.0 ? 1.0 : first);
    for (std::size_t row = column + 2; row < size; ++row)
      reflector[row] = matrix[matrix_offset(row, column, size)];
    Real reflector_norm = 0.0;
    for (std::size_t row = column + 1; row < size; ++row)
      reflector_norm = std::hypot(reflector_norm, reflector[row]);
    for (std::size_t row = column + 1; row < size; ++row)
      reflector[row] /= reflector_norm;

    for (std::size_t j = column; j < size; ++j) {
      Real projection = 0.0;
      for (std::size_t row = column + 1; row < size; ++row)
        projection += reflector[row] * matrix[matrix_offset(row, j, size)];
      for (std::size_t row = column + 1; row < size; ++row)
        matrix[matrix_offset(row, j, size)] -= 2.0 * reflector[row] * projection;
    }
    for (std::size_t row = 0; row < size; ++row) {
      Real projection = 0.0;
      for (std::size_t j = column + 1; j < size; ++j)
        projection += matrix[matrix_offset(row, j, size)] * reflector[j];
      for (std::size_t j = column + 1; j < size; ++j)
        matrix[matrix_offset(row, j, size)] -= 2.0 * projection * reflector[j];
    }
    for (std::size_t row = 0; row < size; ++row) {
      Real projection = 0.0;
      for (std::size_t j = column + 1; j < size; ++j)
        projection += transform[matrix_offset(row, j, size)] * reflector[j];
      for (std::size_t j = column + 1; j < size; ++j)
        transform[matrix_offset(row, j, size)] -= 2.0 * projection * reflector[j];
    }
    for (std::size_t row = column + 2; row < size; ++row)
      matrix[matrix_offset(row, column, size)] = 0.0;
  }
}

[[nodiscard]] RealEigenDecomposition decompose_real_positive_matrix(
    const std::span<const Real> matrix, const std::size_t size) {
  if (matrix.size() != size * size || size == 0)
    throw std::invalid_argument("dry vertical structure matrix shape is invalid");
  std::vector<Real> schur(matrix.begin(), matrix.end());
  auto transform = identity_matrix(size);
  reduce_to_upper_hessenberg(schur, transform, size);
  std::vector<Real> q;
  std::vector<Real> r;
  std::size_t active = size;
  Real matrix_scale = 1.0;
  for (const Real value : matrix)
    matrix_scale = std::max(matrix_scale, std::abs(value));
  constexpr std::size_t maximum_iterations_per_level = 4000;
  std::size_t iterations = 0;
  while (active > 1) {
    const auto subdiagonal = matrix_offset(active - 1, active - 2, size);
    const Real local_scale =
        std::abs(schur[matrix_offset(active - 2, active - 2, size)]) +
        std::abs(schur[matrix_offset(active - 1, active - 1, size)]) + matrix_scale;
    if (std::abs(schur[subdiagonal]) <=
        128.0 * std::numeric_limits<Real>::epsilon() * local_scale) {
      schur[subdiagonal] = 0.0;
      --active;
      iterations = 0;
      continue;
    }
    if (++iterations > maximum_iterations_per_level)
      throw std::runtime_error("dry vertical eigenvalue iteration did not converge");

    const Real a = schur[matrix_offset(active - 2, active - 2, size)];
    const Real b = schur[matrix_offset(active - 2, active - 1, size)];
    const Real c = schur[matrix_offset(active - 1, active - 2, size)];
    const Real d = schur[matrix_offset(active - 1, active - 1, size)];
    const Real discriminant = (a - d) * (a - d) + 4.0 * b * c;
    if (discriminant < 0.0)
      throw std::runtime_error("dry vertical structure has complex eigenvalues");
    const Real root = std::sqrt(discriminant);
    const Real first = 0.5 * (a + d + root);
    const Real second = 0.5 * (a + d - root);
    const Real shift = std::abs(first - d) < std::abs(second - d) ? first : second;

    std::vector<Real> leading(active * active);
    for (std::size_t row = 0; row < active; ++row) {
      for (std::size_t column = 0; column < active; ++column) {
        leading[matrix_offset(row, column, active)] =
            schur[matrix_offset(row, column, size)] - (row == column ? shift : 0.0);
      }
    }
    qr_factorize(leading, active, q, r);
    std::vector<Real> embedded = identity_matrix(size);
    for (std::size_t row = 0; row < active; ++row)
      for (std::size_t column = 0; column < active; ++column)
        embedded[matrix_offset(row, column, size)] =
            q[matrix_offset(row, column, active)];
    std::vector<Real> transposed(size * size);
    for (std::size_t row = 0; row < size; ++row)
      for (std::size_t column = 0; column < size; ++column)
        transposed[matrix_offset(row, column, size)] =
            embedded[matrix_offset(column, row, size)];
    schur =
        multiply_matrices(multiply_matrices(transposed, schur, size), embedded, size);
    transform = multiply_matrices(transform, embedded, size);
  }

  std::vector<Real> eigenvalues(size);
  std::vector<Real> eigenvectors(size * size, 0.0);
  for (std::size_t index = 0; index < size; ++index) {
    const Real eigenvalue = schur[matrix_offset(index, index, size)];
    if (!(eigenvalue > 0.0) || !std::isfinite(eigenvalue))
      throw std::runtime_error("dry vertical structure is not positive wave-like");
    eigenvalues[index] = eigenvalue;
    std::vector<Real> schur_vector(size, 0.0);
    schur_vector[index] = 1.0;
    for (std::size_t reverse = 0; reverse < index; ++reverse) {
      const std::size_t row = index - reverse - 1;
      Real value = 0.0;
      for (std::size_t column = row + 1; column <= index; ++column)
        value += schur[matrix_offset(row, column, size)] * schur_vector[column];
      Real denominator = schur[matrix_offset(row, row, size)] - eigenvalue;
      const Real floor = 128.0 * std::numeric_limits<Real>::epsilon() * matrix_scale;
      if (std::abs(denominator) < floor)
        denominator = std::copysign(floor, denominator == 0.0 ? 1.0 : denominator);
      schur_vector[row] = -value / denominator;
    }
    Real norm_squared = 0.0;
    for (std::size_t row = 0; row < size; ++row) {
      Real value = 0.0;
      for (std::size_t column = 0; column < size; ++column)
        value += transform[matrix_offset(row, column, size)] * schur_vector[column];
      eigenvectors[matrix_offset(row, index, size)] = value;
      norm_squared += value * value;
    }
    const Real vector_norm = std::sqrt(norm_squared);
    if (!(vector_norm > 0.0) || !std::isfinite(vector_norm))
      throw std::runtime_error("dry vertical eigenvector is invalid");
    for (std::size_t row = 0; row < size; ++row)
      eigenvectors[matrix_offset(row, index, size)] /= vector_norm;
  }

  std::vector<std::size_t> order(size);
  std::iota(order.begin(), order.end(), 0);
  std::ranges::sort(order, [&](const std::size_t left, const std::size_t right) {
    return eigenvalues[left] > eigenvalues[right];
  });
  std::vector<Real> sorted_values(size);
  std::vector<Real> sorted_vectors(size * size);
  for (std::size_t mode = 0; mode < size; ++mode) {
    sorted_values[mode] = eigenvalues[order[mode]];
    for (std::size_t level = 0; level < size; ++level)
      sorted_vectors[matrix_offset(level, mode, size)] =
          eigenvectors[matrix_offset(level, order[mode], size)];
  }
  auto inverse = invert_matrix(sorted_vectors, size);
  return {.eigenvalues = std::move(sorted_values),
          .right_eigenvectors = std::move(sorted_vectors),
          .inverse_eigenvectors = std::move(inverse)};
}

}  // namespace

std::span<const Real> DryHydrostaticVerticalModes::eigenvector(
    const std::size_t mode) const {
  if (mode >= mode_count() || eigenvectors.size() != mode_count() * levels)
    throw std::out_of_range("dry vertical mode index or shape is invalid");
  return {eigenvectors.data() + mode * levels, levels};
}

DryHydrostaticReferenceColumn make_dry_hydrostatic_reference_column(
    const AtmosphericHybridCoordinate& coordinate, const PlanetParameters& planet,
    const SemiImplicitParameters& parameters) {
  planet.validate();
  require_positive(parameters.reference_surface_pressure_pa,
                   "semi-implicit reference surface pressure");
  require_positive(parameters.reference_temperature_k,
                   "semi-implicit reference temperature");
  DryHydrostaticReferenceColumn result{
      .surface_pressure_pa = parameters.reference_surface_pressure_pa,
      .temperature_k = parameters.reference_temperature_k,
      .geometry = coordinate.geometry(parameters.reference_surface_pressure_pa,
                                      planet.gravity_m_s2, planet.gas_constant_j_kg_k,
                                      planet.heat_capacity_cp_j_kg_k,
                                      planet.reference_pressure_pa),
      .potential_temperature_k = {},
      .potential_temperature_mass_k_kg_m2 = {},
      .hydrostatic = {}};
  const auto levels = coordinate.levels();
  result.potential_temperature_k.resize(levels);
  result.potential_temperature_mass_k_kg_m2.resize(levels);
  for (std::size_t level = 0; level < levels; ++level) {
    const Real potential_temperature =
        parameters.reference_temperature_k / result.geometry.exner_full[level];
    result.potential_temperature_k[level] = potential_temperature;
    result.potential_temperature_mass_k_kg_m2[level] =
        result.geometry.air_mass_kg_m2[level] * potential_temperature;
  }
  result.hydrostatic = integrate_hydrostatic_column(
      result.geometry, result.potential_temperature_k, planet.heat_capacity_cp_j_kg_k,
      planet.gravity_m_s2, 0.0);
  return result;
}

DryHydrostaticVerticalModes make_dry_hydrostatic_external_mode(
    const DryHydrostaticReferenceColumn& reference, const PlanetParameters& planet) {
  planet.validate();
  const std::size_t levels = reference.geometry.air_mass_kg_m2.size();
  if (levels == 0 || reference.potential_temperature_k.size() != levels ||
      reference.temperature_k <= 0.0 || !std::isfinite(reference.temperature_k))
    throw std::invalid_argument("dry reference column shape or temperature is invalid");
  const Real gamma = planet.heat_capacity_cp_j_kg_k /
                     (planet.heat_capacity_cp_j_kg_k - planet.gas_constant_j_kg_k);
  const Real eigenvalue = gamma * planet.gas_constant_j_kg_k * reference.temperature_k;
  const Real normalization = 1.0 / std::sqrt(static_cast<Real>(levels));
  return {.levels = levels,
          .phase_speed_m_s = {std::sqrt(eigenvalue)},
          .eigenvalue_m2_s2 = {eigenvalue},
          .eigenvectors = std::vector<Real>(levels, normalization),
          .inverse_eigenvectors = std::vector<Real>(levels, normalization),
          .vertical_structure_m2_s2 = {}};
}

DryHydrostaticVerticalModes make_dry_hydrostatic_vertical_modes(
    const DryHydrostaticReferenceColumn& reference, const PlanetParameters& planet,
    const DryHydrostaticFastOperator& op) {
  planet.validate();
  const std::size_t levels = op.levels;
  if (levels == 0 || reference.geometry.air_mass_kg_m2.size() != levels ||
      op.b_half.size() != levels + 1 || op.reference_air_mass_kg_m2.size() != levels ||
      op.reference_potential_temperature_k.size() != levels ||
      op.reference_interface_potential_temperature_k.size() != levels + 1 ||
      op.reference_specific_volume_m3_kg.size() != levels ||
      op.pressure_from_surface_pressure.size() != levels ||
      op.geopotential_from_surface_pressure.size() != levels ||
      op.geopotential_from_potential_temperature_mass.size() != levels * levels)
    throw std::invalid_argument("dry vertical structure inputs are malformed");

  std::vector<Real> structure(levels * levels, 0.0);
  std::vector<Real> horizontal_mass(levels, 0.0);
  VerticalMassFlux mass_flux;
  std::vector<Real> thermal_tendency(levels, 0.0);
  for (std::size_t source = 0; source < levels; ++source) {
    std::fill(horizontal_mass.begin(), horizontal_mass.end(), 0.0);
    horizontal_mass[source] = -1.0;
    diagnose_vertical_mass_flux(horizontal_mass, op.b_half, planet.gravity_m_s2,
                                mass_flux);
    for (std::size_t level = 0; level < levels; ++level) {
      thermal_tendency[level] =
          op.reference_potential_temperature_k[level] * horizontal_mass[level] +
          mass_flux.interface_flux_kg_m2_s[level] *
              op.reference_interface_potential_temperature_k[level] -
          mass_flux.interface_flux_kg_m2_s[level + 1] *
              op.reference_interface_potential_temperature_k[level + 1];
    }
    for (std::size_t level = 0; level < levels; ++level) {
      Real pressure_potential = (op.geopotential_from_surface_pressure[level] +
                                 op.reference_specific_volume_m3_kg[level] *
                                     op.pressure_from_surface_pressure[level]) *
                                mass_flux.surface_pressure_tendency_pa_s;
      for (std::size_t thermal_level = 0; thermal_level < levels; ++thermal_level) {
        pressure_potential +=
            op.geopotential_from_potential_temperature_mass[matrix_offset(
                level, thermal_level, levels)] *
            thermal_tendency[thermal_level];
      }
      structure[matrix_offset(level, source, levels)] =
          -op.reference_air_mass_kg_m2[level] * pressure_potential;
    }
  }

  const auto decomposition = decompose_real_positive_matrix(structure, levels);
  DryHydrostaticVerticalModes result{
      .levels = levels,
      .phase_speed_m_s = std::vector<Real>(levels),
      .eigenvalue_m2_s2 = decomposition.eigenvalues,
      .eigenvectors = std::vector<Real>(levels * levels),
      .inverse_eigenvectors = decomposition.inverse_eigenvectors,
      .vertical_structure_m2_s2 = std::move(structure)};
  // Convert column-major right eigenvectors to the public mode-major layout.
  for (std::size_t mode = 0; mode < levels; ++mode) {
    result.phase_speed_m_s[mode] = std::sqrt(result.eigenvalue_m2_s2[mode]);
    for (std::size_t level = 0; level < levels; ++level)
      result.eigenvectors[matrix_offset(mode, level, levels)] =
          decomposition.right_eigenvectors[matrix_offset(level, mode, levels)];
  }
  return result;
}

std::vector<std::size_t> select_implicit_vertical_modes(
    const CubedSphereGrid& grid, const DryHydrostaticVerticalModes& modes,
    const Real time_step_s, const Real wave_cfl_threshold,
    const std::size_t maximum_implicit_modes) {
  require_positive(time_step_s, "semi-implicit mode-selection time step");
  require_positive(wave_cfl_threshold, "semi-implicit wave CFL threshold");
  if (maximum_implicit_modes == 0)
    throw std::invalid_argument("maximum implicit mode count must be positive");
  if (modes.mode_count() == 0 || modes.eigenvalue_m2_s2.size() != modes.mode_count() ||
      modes.eigenvectors.size() != modes.mode_count() * modes.levels ||
      (!modes.inverse_eigenvectors.empty() &&
       modes.inverse_eigenvectors.size() != modes.mode_count() * modes.levels))
    throw std::invalid_argument("dry vertical modes are empty or malformed");
  Real maximum_geometry_factor_m_inverse = 0.0;
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
    Real perimeter = 0.0;
    for (const auto& edge : grid.cell_cache()[cell].edges)
      perimeter += grid.edges()[edge.edge].length_m;
    maximum_geometry_factor_m_inverse = std::max(
        maximum_geometry_factor_m_inverse, perimeter / grid.cells()[cell].area_m2);
  }
  std::vector<std::size_t> selected;
  Real previous_speed = std::numeric_limits<Real>::infinity();
  for (std::size_t mode = 0; mode < modes.mode_count(); ++mode) {
    const Real speed = modes.phase_speed_m_s[mode];
    if (!(speed > 0.0) || !std::isfinite(speed) || speed > previous_speed)
      throw std::invalid_argument(
          "dry vertical mode speeds must be positive and descending");
    previous_speed = speed;
    const Real courant = time_step_s * speed * maximum_geometry_factor_m_inverse;
    if (courant > wave_cfl_threshold) selected.push_back(mode);
  }
  if (selected.size() > maximum_implicit_modes)
    throw std::runtime_error(
        "required implicit vertical modes exceed the configured maximum");
  return selected;
}

Real maximum_vertical_mode_courant(const CubedSphereGrid& grid,
                                   const DryHydrostaticVerticalModes& modes,
                                   const Real time_step_s) {
  require_positive(time_step_s, "semi-implicit wave Courant time step");
  if (modes.mode_count() == 0 || !(modes.phase_speed_m_s.front() > 0.0) ||
      !std::isfinite(modes.phase_speed_m_s.front()))
    throw std::invalid_argument("dry vertical modes are empty or malformed");
  Real maximum_geometry_factor_m_inverse = 0.0;
  for (std::size_t cell = 0; cell < grid.cell_count(); ++cell) {
    Real perimeter = 0.0;
    for (const auto& edge : grid.cell_cache()[cell].edges)
      perimeter += grid.edges()[edge.edge].length_m;
    maximum_geometry_factor_m_inverse = std::max(
        maximum_geometry_factor_m_inverse, perimeter / grid.cells()[cell].area_m2);
  }
  return time_step_s * modes.phase_speed_m_s.front() *
         maximum_geometry_factor_m_inverse;
}

void write_dry_hydrostatic_vertical_mode_metadata(
    std::ostream& output, const DryHydrostaticVerticalModes& modes) {
  if (modes.mode_count() == 0 || modes.eigenvalue_m2_s2.size() != modes.mode_count())
    throw std::invalid_argument("dry vertical modes are empty or malformed");
  output << std::setprecision(std::numeric_limits<Real>::max_digits10)
         << "semi_implicit.vertical_mode_count = " << modes.mode_count() << '\n';
  for (std::size_t mode = 0; mode < modes.mode_count(); ++mode) {
    output << "semi_implicit.mode_" << mode
           << "_phase_speed_m_s = " << modes.phase_speed_m_s[mode] << '\n';
  }
  if (!output)
    throw std::runtime_error("failed while writing dry vertical mode metadata");
}

void project_onto_vertical_modes(const DryHydrostaticVerticalModes& modes,
                                 const std::span<const Real> level_values,
                                 const std::span<Real> mode_values) {
  if (level_values.size() != modes.levels || mode_values.size() != modes.mode_count())
    throw std::invalid_argument("vertical mode projection shapes differ");
  for (std::size_t mode = 0; mode < modes.mode_count(); ++mode) {
    mode_values[mode] = 0.0;
    const auto projection =
        modes.inverse_eigenvectors.empty()
            ? modes.eigenvector(mode)
            : std::span<const Real>(
                  modes.inverse_eigenvectors.data() + mode * modes.levels,
                  modes.levels);
    for (std::size_t level = 0; level < modes.levels; ++level)
      mode_values[mode] += projection[level] * level_values[level];
  }
}

void reconstruct_from_vertical_modes(const DryHydrostaticVerticalModes& modes,
                                     const std::span<const Real> mode_values,
                                     const std::span<Real> level_values) {
  if (level_values.size() != modes.levels || mode_values.size() != modes.mode_count())
    throw std::invalid_argument("vertical mode reconstruction shapes differ");
  std::fill(level_values.begin(), level_values.end(), 0.0);
  for (std::size_t mode = 0; mode < modes.mode_count(); ++mode) {
    const auto eigenvector = modes.eigenvector(mode);
    for (std::size_t level = 0; level < modes.levels; ++level)
      level_values[level] += mode_values[mode] * eigenvector[level];
  }
}

}  // namespace mps
