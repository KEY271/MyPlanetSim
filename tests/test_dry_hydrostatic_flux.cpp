#include "myplanetsim/dynamics/dry_hydrostatic_flux.hpp"
#include "support/test.hpp"
MPS_TEST_CASE("layer flux preserves constant mass-specific scalars"){
 mps::EdgeTangentBasis b{.normal={1,0,0},.tangent={0,1,0}}; mps::DryHydrostaticPrimitive l{2,{3,1,0},300,0.25,280},r{4,{1,-1,0},300,0.25,290}; auto f=mps::rusanov_dry_hydrostatic_flux(l,r,b,287,1004);
 MPS_CHECK_NEAR(f.potential_temperature_mass_k_kg_m_s,300*f.air_mass_kg_m_s,1e-10); MPS_CHECK_NEAR(f.tracer_mass_kg_m_s,0.25*f.air_mass_kg_m_s,1e-12); MPS_CHECK(f.maximum_wave_speed_m_s>0);
}
int main(){return mps::test::run_all();}
