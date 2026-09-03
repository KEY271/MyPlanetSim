#pragma once

#include <cstddef>
#include <cstdint>
#include <limits>

namespace mps {

using Real = double;
using Index = std::ptrdiff_t;
using Seed = std::uint64_t;

static_assert(sizeof(Real) == 8, "MyPlanetSim requires 64-bit double precision");
static_assert(std::numeric_limits<Real>::is_iec559,
              "MyPlanetSim requires IEEE 754 floating-point semantics");
static_assert(std::numeric_limits<Real>::digits == 53,
              "MyPlanetSim requires IEEE 754 binary64 precision");

}  // namespace mps
