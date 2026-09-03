#include <cmath>
#include <stdexcept>

#include "myplanetsim/dynamics/shallow_water_flux.hpp"
#include "support/test.hpp"

namespace {

constexpr mps::EdgeTangentBasis kBasis{.normal = {1.0, 0.0, 0.0},
                                       .tangent = {0.0, 1.0, 0.0}};

}  // namespace

MPS_TEST_CASE("Rusanov flux recovers the physical flux for identical states") {
  constexpr mps::ShallowWaterPrimitive state{.depth_m = 4.0,
                                             .velocity_m_s = {3.0, 2.0, 0.0}};
  const auto flux = mps::rusanov_shallow_water_flux(state, state, kBasis, 5.0);
  MPS_CHECK_NEAR(flux.mass_m2_s, 12.0, 0.0);
  MPS_CHECK_NEAR(flux.momentum_m3_s2.x, 76.0, 0.0);
  MPS_CHECK_NEAR(flux.momentum_m3_s2.y, 24.0, 0.0);
  MPS_CHECK_NEAR(flux.momentum_m3_s2.z, 0.0, 0.0);
  MPS_CHECK_NEAR(flux.maximum_wave_speed_m_s, 3.0 + std::sqrt(20.0), 0.0);
}

MPS_TEST_CASE("Rusanov lake-at-rest flux contains pressure but no mass flux") {
  constexpr mps::ShallowWaterPrimitive state{.depth_m = 2.0,
                                             .velocity_m_s = {0.0, 0.0, 0.0}};
  const auto flux = mps::rusanov_shallow_water_flux(state, state, kBasis, 3.0);
  MPS_CHECK_EQ(flux.mass_m2_s, 0.0);
  MPS_CHECK_NEAR(flux.momentum_m3_s2.x, 6.0, 0.0);
  MPS_CHECK_NEAR(flux.momentum_m3_s2.y, 0.0, 0.0);
}

MPS_TEST_CASE("Rusanov flux changes sign when edge orientation is reversed") {
  constexpr mps::ShallowWaterPrimitive left{.depth_m = 2.0,
                                            .velocity_m_s = {1.0, 0.5, 0.0}};
  constexpr mps::ShallowWaterPrimitive right{.depth_m = 3.0,
                                             .velocity_m_s = {-0.25, 2.0, 0.0}};
  const auto forward = mps::rusanov_shallow_water_flux(left, right, kBasis, 4.0);
  constexpr mps::EdgeTangentBasis reversed{.normal = {-1.0, 0.0, 0.0},
                                           .tangent = {0.0, -1.0, 0.0}};
  const auto backward = mps::rusanov_shallow_water_flux(right, left, reversed, 4.0);
  MPS_CHECK_NEAR(backward.mass_m2_s, -forward.mass_m2_s, 2.0e-15);
  MPS_CHECK_NEAR(mps::norm(backward.momentum_m3_s2 + forward.momentum_m3_s2), 0.0,
                 2.0e-14);
}

MPS_TEST_CASE("Rusanov flux rejects nonphysical edge states") {
  constexpr mps::ShallowWaterPrimitive invalid{.depth_m = 0.0, .velocity_m_s = {}};
  MPS_CHECK_THROWS_AS(mps::rusanov_shallow_water_flux(invalid, invalid, kBasis, 1.0),
                      std::invalid_argument);
}

int main() { return mps::test::run_all(); }
