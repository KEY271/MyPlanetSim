#include "myplanetsim/dynamics/dry_hydrostatic_benchmarks.hpp"
#include "support/test.hpp"
MPS_TEST_CASE("phase 5 benchmark presets produce shaped finite states") {
  mps::ExperimentConfig c{.kind=mps::ExperimentKind::kDryHydrostatic,.planet={2,0,10,287,1004,100000},.run={0,1,.1,0},.grid={2},.vertical={.levels=2,.a_half_pa={1000,500,0},.b_half={0,.5,1},.surface_pressure_pa=100000,.minimum_surface_pressure_pa=90000,.maximum_surface_pressure_pa=110000,.minimum_pressure_thickness_pa=100,.initial_temperature_k=280,.initial_potential_temperature_k=300,.temperature_floor_k=100,.transport_scheme=mps::VerticalTransportScheme::kDonorCell,.limiter=mps::VerticalLimiterKind::kNone,.cfl=.5},.dry_hydrostatic={},.diagnostics={1},.output_directory="x"};
  mps::CubedSphereGrid grid(2,2);mps::AtmosphericHybridCoordinate z({c.vertical.a_half_pa,c.vertical.b_half},90000,110000,100);
  for(auto kind:{mps::DryHydrostaticTestCase::kIsothermalRest,mps::DryHydrostaticTestCase::kSolidBodyTransport,mps::DryHydrostaticTestCase::kDcmipDeformational,mps::DryHydrostaticTestCase::kDcmipHadley,mps::DryHydrostaticTestCase::kLinearWave,mps::DryHydrostaticTestCase::kUmjs14Steady,mps::DryHydrostaticTestCase::kUmjs14Baroclinic}){c.dry_hydrostatic.test_case=kind;auto s=mps::initialize_dry_hydrostatic_benchmark(c,grid,z);MPS_CHECK_EQ(s.surface_pressure_pa.size(),grid.cell_count());for(auto v:mps::flatten_dry_hydrostatic_state(s,2))MPS_CHECK(std::isfinite(v));}
}
int main(){return mps::test::run_all();}
