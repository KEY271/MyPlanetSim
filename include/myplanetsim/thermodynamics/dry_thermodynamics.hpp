#pragma once

#include "myplanetsim/core/types.hpp"

namespace mps::thermodynamics {

[[nodiscard]] Real exner(Real pressure_pa, Real gas_constant_j_kg_k,
                         Real heat_capacity_cp_j_kg_k, Real reference_pressure_pa);
[[nodiscard]] Real temperature_from_potential_temperature(Real potential_temperature_k,
                                                          Real pressure_pa,
                                                          Real gas_constant_j_kg_k,
                                                          Real heat_capacity_cp_j_kg_k,
                                                          Real reference_pressure_pa);
[[nodiscard]] Real potential_temperature_from_temperature(Real temperature_k,
                                                          Real pressure_pa,
                                                          Real gas_constant_j_kg_k,
                                                          Real heat_capacity_cp_j_kg_k,
                                                          Real reference_pressure_pa);
[[nodiscard]] Real density(Real pressure_pa, Real temperature_k,
                           Real gas_constant_j_kg_k);

}  // namespace mps::thermodynamics
