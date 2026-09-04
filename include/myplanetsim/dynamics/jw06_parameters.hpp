#pragma once

#include "myplanetsim/core/planet_parameters.hpp"

namespace mps {

inline constexpr PlanetParameters kJw06Planet{6371229.0, 7.29212e-5, 9.80616,
                                              287.0,     1004.5,     100000.0};
inline constexpr Real kJw06U0Mps = 35.0;
inline constexpr Real kJw06Eta0 = 0.252;

}  // namespace mps
