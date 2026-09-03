#include "myplanetsim/transport/spherical_transport.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>

#include "myplanetsim/numerics/explicit_steppers.hpp"
#include "myplanetsim/numerics/spherical_operators.hpp"

namespace mps {
namespace {

[[nodiscard]] Vec3 rotate_about_axis(const Vec3 point, const Vec3 unit_axis,
                                     const Real angle) noexcept {
  return std::cos(angle) * point + std::sin(angle) * cross(unit_axis, point) +
         (1.0 - std::cos(angle)) * dot(unit_axis, point) * unit_axis;
}

[[nodiscard]] Vec3 face_displacement(const Vec3 center, const Vec3 face,
                                     const Real radius_m) {
  const Real cosine = std::clamp(dot(center, face), -1.0, 1.0);
  const Real angle = std::acos(cosine);
  const Real sine = std::sin(angle);
  if (!(sine > 0.0)) {
    return {};
  }
  return (radius_m * angle / sine) * (face - cosine * center);
}

[[nodiscard]] Real streamfunction(const TransportTestCase test_case, const Vec3 point,
                                  const Real time_s, const Real period_s,
                                  const Vec3 axis, const Real angular_speed_rad_s,
                                  const Real radius_m) {
  const Real solid = radius_m * radius_m * angular_speed_rad_s * dot(axis, point);
  if (test_case == TransportTestCase::kSolidBody) {
    return solid;
  }
  const Real phase =
      period_s > 0.0 ? std::cos(std::numbers::pi_v<Real> * time_s / period_s) : 1.0;
  const Real longitude = std::atan2(point.y, point.x);
  const Real deformation = 0.25 * radius_m * radius_m * angular_speed_rad_s *
                           std::sin(2.0 * longitude) * (1.0 - point.z * point.z) *
                           phase;
  return solid + deformation;
}

void validate_tracer(const CubedSphereGrid& grid, const std::span<const Real> tracer) {
  if (tracer.size() != grid.cell_count()) {
    throw std::invalid_argument("tracer size does not match cubed-sphere grid");
  }
  if (!std::ranges::all_of(tracer,
                           [](const Real value) { return std::isfinite(value); })) {
    throw std::invalid_argument("tracer contains a non-finite value");
  }
}

[[nodiscard]] TransportDiagnostics diagnose(const CubedSphereGrid& grid,
                                            const std::span<const Real> initial,
                                            const std::span<const Real> final,
                                            const std::span<const Real> exact,
                                            const std::uint64_t limiter_activations) {
  const Real initial_mass = tracer_mass(grid, initial);
  const Real final_mass = tracer_mass(grid, final);
  Real weighted_absolute = 0.0;
  Real weighted_squared = 0.0;
  Real linf = 0.0;
  Real area = 0.0;
  Real seam_error = 0.0;
  Real seam_area = 0.0;
  Real corner_error = 0.0;
  Real corner_area = 0.0;
  const auto [initial_minimum, initial_maximum] =
      std::minmax_element(initial.begin(), initial.end());
  Real below = 0.0;
  Real above = 0.0;
  Real interior_mixing = 0.0;
  for (std::size_t index = 0; index < final.size(); ++index) {
    const Real error = final[index] - exact[index];
    const Real weight = grid.cells()[index].area_m2;
    weighted_absolute += weight * std::abs(error);
    weighted_squared += weight * error * error;
    linf = std::max(linf, std::abs(error));
    area += weight;
    const auto id = grid.cells()[index].id;
    const bool i_edge = id.i == 0 || id.i == grid.cells_per_panel() - 1;
    const bool j_edge = id.j == 0 || id.j == grid.cells_per_panel() - 1;
    if (i_edge || j_edge) {
      seam_error += weight * error * error;
      seam_area += weight;
    }
    if (i_edge && j_edge) {
      corner_error += weight * error * error;
      corner_area += weight;
    }
    below += weight * std::max(*initial_minimum - final[index], 0.0);
    above += weight * std::max(final[index] - *initial_maximum, 0.0);
    if (final[index] >= *initial_minimum && final[index] <= *initial_maximum) {
      interior_mixing += weight * std::abs(error);
    }
  }
  const auto [minimum, maximum] = std::minmax_element(final.begin(), final.end());
  Real filament = 0.0;
  for (int level = 1; level <= 9; ++level) {
    const Real threshold = static_cast<Real>(level) / 10.0;
    Real initial_area = 0.0;
    Real final_area = 0.0;
    for (std::size_t index = 0; index < final.size(); ++index) {
      const Real weight = grid.cells()[index].area_m2;
      initial_area += initial[index] >= threshold ? weight : 0.0;
      final_area += final[index] >= threshold ? weight : 0.0;
    }
    if (initial_area > 0.0) {
      filament += std::min(final_area / initial_area, 1.0);
    } else {
      filament += 1.0;
    }
  }
  filament /= 9.0;
  return {initial_mass,
          final_mass,
          (final_mass - initial_mass) /
              std::max(std::abs(initial_mass), std::numeric_limits<Real>::min()),
          weighted_absolute / area,
          std::sqrt(weighted_squared / area),
          linf,
          *minimum,
          *maximum,
          std::sqrt(seam_error / seam_area),
          std::sqrt(corner_error / corner_area),
          filament,
          below / area,
          above / area,
          interior_mixing / area,
          limiter_activations};
}

}  // namespace

Vec3 solid_body_velocity(const Vec3 unit_position, const Real radius_m,
                         const Vec3 rotation_axis, const Real angular_speed_rad_s) {
  if (!(radius_m > 0.0) || !std::isfinite(radius_m) ||
      !std::isfinite(angular_speed_rad_s)) {
    throw std::invalid_argument("solid-body wind parameters must be finite and valid");
  }
  const Vec3 position = normalize(unit_position);
  const Vec3 axis = normalize(rotation_axis);
  return radius_m * angular_speed_rad_s * cross(axis, position);
}

std::vector<Real> prescribed_edge_fluxes(const CubedSphereGrid& grid,
                                         const TransportTestCase test_case,
                                         const Real time_s, const Real period_s,
                                         const Vec3 rotation_axis,
                                         const Real angular_speed_rad_s) {
  const Vec3 axis = normalize(rotation_axis);
  std::vector<Real> fluxes(grid.edge_count());
  for (const auto& edge : grid.edges()) {
    const Vec3 first = grid.vertices()[edge.first_vertex].position;
    const Vec3 second = grid.vertices()[edge.second_vertex].position;
    const Real nondivergent = streamfunction(test_case, second, time_s, period_s, axis,
                                             angular_speed_rad_s, grid.radius_m()) -
                              streamfunction(test_case, first, time_s, period_s, axis,
                                             angular_speed_rad_s, grid.radius_m());
    if (test_case != TransportTestCase::kDivergent) {
      fluxes[edge.id] = nondivergent;
    } else {
      const Real phase =
          period_s > 0.0 ? std::cos(std::numbers::pi_v<Real> * time_s / period_s) : 1.0;
      const Vec3 potential_velocity = grid.radius_m() * angular_speed_rad_s * phase *
                                      project_tangent(Vec3{0.0, 0.0, 1.0}, edge.center);
      fluxes[edge.id] =
          nondivergent +
          dot(potential_velocity, edge.outward_normal_from_left) * edge.length_m;
    }
  }
  return fluxes;
}

Real initial_tracer_value(const InitialConditionKind kind, const Vec3 unit_position) {
  const Vec3 point = normalize(unit_position);
  constexpr Vec3 center{1.0, 0.0, 0.0};
  const Real distance = safe_angle(center, point);
  switch (kind) {
    case InitialConditionKind::kConstant:
      return 1.0;
    case InitialConditionKind::kGaussianHill:
      return std::exp(-(1.0 - dot(center, point)));
    case InitialConditionKind::kCosineBell: {
      constexpr Real radius = std::numbers::pi_v<Real> / 3.0;
      return distance < radius
                 ? 0.5 * (1.0 + std::cos(std::numbers::pi_v<Real> * distance / radius))
                 : 0.0;
    }
    case InitialConditionKind::kSlottedCylinder: {
      constexpr Real radius = 0.5;
      if (distance >= radius) {
        return 0.0;
      }
      const Real longitude = std::atan2(point.y, point.x);
      return std::abs(longitude) < 0.08 && point.z < 0.25 ? 0.0 : 1.0;
    }
  }
  throw std::logic_error("unknown initial condition");
}

std::vector<Real> make_initial_tracer(const CubedSphereGrid& grid,
                                      const InitialConditionKind kind) {
  std::vector<Real> tracer(grid.cell_count());
  for (std::size_t index = 0; index < tracer.size(); ++index) {
    tracer[index] = initial_tracer_value(kind, grid.cells()[index].center);
  }
  return tracer;
}

std::vector<Real> exact_solid_body_tracer(const CubedSphereGrid& grid,
                                          const InitialConditionKind kind,
                                          const Vec3 rotation_axis,
                                          const Real angular_speed_rad_s,
                                          const Real elapsed_time_s) {
  const Vec3 axis = normalize(rotation_axis);
  std::vector<Real> tracer(grid.cell_count());
  for (std::size_t index = 0; index < tracer.size(); ++index) {
    const Vec3 departure = rotate_about_axis(grid.cells()[index].center, axis,
                                             -angular_speed_rad_s * elapsed_time_s);
    tracer[index] = initial_tracer_value(kind, departure);
  }
  return tracer;
}

Real tracer_mass(const CubedSphereGrid& grid, const std::span<const Real> tracer) {
  validate_tracer(grid, tracer);
  Real mass = 0.0;
  for (std::size_t index = 0; index < tracer.size(); ++index) {
    mass += grid.cells()[index].area_m2 * tracer[index];
  }
  return mass;
}

Real stable_transport_time_step(const CubedSphereGrid& grid,
                                const std::span<const Real> oriented_edge_flux,
                                const Real cfl, const Real maximum_time_step_s) {
  if (oriented_edge_flux.size() != grid.edge_count() || !(cfl > 0.0 && cfl <= 1.0) ||
      !(maximum_time_step_s > 0.0) || !std::isfinite(maximum_time_step_s)) {
    throw std::invalid_argument("invalid transport CFL arguments");
  }
  std::vector<Real> outgoing(grid.cell_count(), 0.0);
  for (const auto& edge : grid.edges()) {
    const Real flux = oriented_edge_flux[edge.id];
    if (!std::isfinite(flux)) {
      throw std::invalid_argument("edge flux contains a non-finite value");
    }
    outgoing[grid.cell_index(edge.left_cell)] += std::max(flux, 0.0);
    outgoing[grid.cell_index(edge.right_cell)] += std::max(-flux, 0.0);
  }
  Real time_step = maximum_time_step_s;
  for (std::size_t index = 0; index < outgoing.size(); ++index) {
    if (outgoing[index] > 0.0) {
      time_step =
          std::min(time_step, cfl * grid.cells()[index].area_m2 / outgoing[index]);
    }
  }
  return time_step;
}

SphericalTransportRhs::SphericalTransportRhs(const CubedSphereGrid& grid,
                                             const TransportParameters& parameters,
                                             const Real period_s)
    : grid_(grid),
      parameters_(parameters),
      period_s_(period_s),
      axis_(normalize({parameters.rotation_axis_x, parameters.rotation_axis_y,
                       parameters.rotation_axis_z})) {}

void SphericalTransportRhs::operator()(const Real time_s,
                                       const std::span<const Real> tracer,
                                       const std::span<Real> tendency) {
  validate_tracer(grid_, tracer);
  if (tendency.size() != grid_.cell_count()) {
    throw std::invalid_argument("transport tendency size does not match grid");
  }
  const auto fluxes =
      prescribed_edge_fluxes(grid_, parameters_.test_case, time_s, period_s_, axis_,
                             parameters_.angular_speed_rad_s);
  std::vector<Vec3> gradients;
  std::vector<Real> factors(grid_.cell_count(), 1.0);
  if (parameters_.scheme == TransportScheme::kLinear) {
    gradients = least_squares_gradient(grid_, tracer);
    if (parameters_.limiter == LimiterKind::kBarthJespersen) {
      for (std::size_t index = 0; index < grid_.cell_count(); ++index) {
        const auto& cell = grid_.cells()[index];
        Real minimum = tracer[index];
        Real maximum = tracer[index];
        for (const auto edge_id : grid_.cell_edges(cell.id)) {
          const Real neighbor =
              tracer[grid_.cell_index(grid_.neighbor_across(edge_id, cell.id))];
          minimum = std::min(minimum, neighbor);
          maximum = std::max(maximum, neighbor);
        }
        for (const auto edge_id : grid_.cell_edges(cell.id)) {
          const Real increment =
              dot(gradients[index],
                  face_displacement(cell.center, grid_.edge(edge_id).center,
                                    grid_.radius_m()));
          if (increment > 0.0) {
            factors[index] =
                std::min(factors[index], (maximum - tracer[index]) / increment);
          } else if (increment < 0.0) {
            factors[index] =
                std::min(factors[index], (minimum - tracer[index]) / increment);
          }
        }
        factors[index] = std::clamp(factors[index], 0.0, 1.0);
        if (factors[index] < 1.0 - 1.0e-14) {
          ++limiter_activations_;
        }
      }
    }
  }

  std::vector<Real> tracer_flux(grid_.edge_count());
  for (const auto& edge : grid_.edges()) {
    const auto left = grid_.cell_index(edge.left_cell);
    const auto right = grid_.cell_index(edge.right_cell);
    Real left_value = tracer[left];
    Real right_value = tracer[right];
    if (parameters_.scheme == TransportScheme::kLinear) {
      left_value +=
          factors[left] *
          dot(gradients[left], face_displacement(grid_.cells()[left].center,
                                                 edge.center, grid_.radius_m()));
      right_value +=
          factors[right] *
          dot(gradients[right], face_displacement(grid_.cells()[right].center,
                                                  edge.center, grid_.radius_m()));
    }
    tracer_flux[edge.id] = fluxes[edge.id] >= 0.0 ? fluxes[edge.id] * left_value
                                                  : fluxes[edge.id] * right_value;
  }
  scatter_oriented_edge_flux(grid_, tracer_flux, tendency);
  for (std::size_t index = 0; index < tendency.size(); ++index) {
    tendency[index] /= grid_.cells()[index].area_m2;
  }
}

TransportResult run_spherical_transport(
    const ExperimentConfig& config, std::optional<TransportState> initial_state,
    const std::optional<std::uint64_t> stop_after_step) {
  config.validate();
  if (config.kind != ExperimentKind::kSphereTransport) {
    throw std::invalid_argument("sphere transport requires sphere_transport config");
  }
  CubedSphereGrid grid(config.grid.cells_per_panel, config.planet.radius_m);
  const auto reference_initial =
      make_initial_tracer(grid, config.transport.initial_condition);
  TransportState state =
      initial_state.has_value()
          ? std::move(*initial_state)
          : TransportState{config.run.start_time_s, 0, reference_initial};
  validate_tracer(grid, state.tracer);
  if (state.time_s < config.run.start_time_s || state.time_s > config.run.end_time_s) {
    throw std::invalid_argument("transport restart time is outside configured run");
  }
  const Real period = config.run.end_time_s - config.run.start_time_s;
  SphericalTransportRhs rhs(grid, config.transport, period);
  SspRk3 stepper(grid.cell_count());
  while (state.time_s < config.run.end_time_s &&
         (!stop_after_step.has_value() || state.step < *stop_after_step)) {
    const auto fluxes = prescribed_edge_fluxes(
        grid, config.transport.test_case, state.time_s, period,
        {config.transport.rotation_axis_x, config.transport.rotation_axis_y,
         config.transport.rotation_axis_z},
        config.transport.angular_speed_rad_s);
    Real time_step = stable_transport_time_step(grid, fluxes, config.transport.cfl,
                                                config.run.time_step_s);
    time_step = std::min(time_step, config.run.end_time_s - state.time_s);
    stepper.step(state.time_s, time_step, state.tracer, rhs);
    state.time_s += time_step;
    ++state.step;
  }
  const bool reached_end = state.time_s == config.run.end_time_s;
  std::vector<Real> exact = reference_initial;
  if (config.transport.test_case == TransportTestCase::kSolidBody) {
    exact = exact_solid_body_tracer(
        grid, config.transport.initial_condition,
        {config.transport.rotation_axis_x, config.transport.rotation_axis_y,
         config.transport.rotation_axis_z},
        config.transport.angular_speed_rad_s, state.time_s - config.run.start_time_s);
  }
  const auto diagnostics =
      diagnose(grid, reference_initial, state.tracer, exact, rhs.limiter_activations());
  return {std::move(state), diagnostics, reached_end};
}

}  // namespace mps
