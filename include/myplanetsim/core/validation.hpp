#pragma once

#include <string_view>

#include "myplanetsim/core/types.hpp"

namespace mps {

void require_finite(Real value, std::string_view parameter_name);
void require_positive(Real value, std::string_view parameter_name);
void require_non_negative(Real value, std::string_view parameter_name);
void require_in_closed_range(Real value, Real lower, Real upper,
                             std::string_view parameter_name);

}  // namespace mps
