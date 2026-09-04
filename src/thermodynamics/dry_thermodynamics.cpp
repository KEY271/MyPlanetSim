#include "myplanetsim/thermodynamics/dry_thermodynamics.hpp"

#include <cmath>
#include <stdexcept>

#include "myplanetsim/core/validation.hpp"

namespace mps::thermodynamics {
namespace {
void validate_constants(const Real gas_constant_j_kg_k,
                        const Real heat_capacity_cp_j_kg_k,
                        const Real reference_pressure_pa) {
  require_positive(gas_constant_j_kg_k, "gas constant");
  require_positive(heat_capacity_cp_j_kg_k, "heat capacity");
  require_positive(reference_pressure_pa, "reference pressure");
  if (!(heat_capacity_cp_j_kg_k > gas_constant_j_kg_k)) {
    throw std::invalid_argument("heat capacity must exceed gas constant");
  }
}
}  // namespace

Real exner(const Real pressure_pa, const Real gas_constant_j_kg_k,
           const Real heat_capacity_cp_j_kg_k, const Real reference_pressure_pa) {
  require_positive(pressure_pa, "pressure");
  validate_constants(gas_constant_j_kg_k, heat_capacity_cp_j_kg_k,
                     reference_pressure_pa);
  const Real value = std::pow(pressure_pa / reference_pressure_pa,
                              gas_constant_j_kg_k / heat_capacity_cp_j_kg_k);
  if (!(value > 0.0) || !std::isfinite(value)) {
    throw std::overflow_error("Exner function is not finite and positive");
  }
  return value;
}

Real temperature_from_potential_temperature(const Real potential_temperature_k,
                                            const Real pressure_pa,
                                            const Real gas_constant_j_kg_k,
                                            const Real heat_capacity_cp_j_kg_k,
                                            const Real reference_pressure_pa) {
  require_positive(potential_temperature_k, "potential temperature");
  return potential_temperature_k * exner(pressure_pa, gas_constant_j_kg_k,
                                         heat_capacity_cp_j_kg_k,
                                         reference_pressure_pa);
}

Real potential_temperature_from_temperature(const Real temperature_k,
                                            const Real pressure_pa,
                                            const Real gas_constant_j_kg_k,
                                            const Real heat_capacity_cp_j_kg_k,
                                            const Real reference_pressure_pa) {
  require_positive(temperature_k, "temperature");
  return temperature_k / exner(pressure_pa, gas_constant_j_kg_k,
                               heat_capacity_cp_j_kg_k, reference_pressure_pa);
}

Real density(const Real pressure_pa, const Real temperature_k,
             const Real gas_constant_j_kg_k) {
  require_positive(pressure_pa, "pressure");
  require_positive(temperature_k, "temperature");
  require_positive(gas_constant_j_kg_k, "gas constant");
  return pressure_pa / (gas_constant_j_kg_k * temperature_k);
}

}  // namespace mps::thermodynamics
