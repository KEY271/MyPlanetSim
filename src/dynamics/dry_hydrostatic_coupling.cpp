#include "myplanetsim/dynamics/dry_hydrostatic_coupling.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
namespace mps {
DryHydrostaticCoupling couple_dry_hydrostatic_columns(const DryHydrostaticState& s,const DryHydrostaticDerived& d,const DryHydrostaticTransportTendency& h,const std::span<const Real> b,const Real g,const VerticalTransportScheme scheme,const VerticalLimiterKind limiter){
 const auto n=d.cells*d.levels;if(h.air_mass.size()!=n||h.momentum.size()!=n||h.potential_temperature_mass.size()!=n||h.tracer_mass.size()!=n||b.size()!=d.levels+1)throw std::invalid_argument("dry coupling shape mismatch"); DryHydrostaticCoupling out;out.surface_pressure_pa_s.resize(d.cells);out.interface_mass_flux_kg_m2_s.resize(d.cells*(d.levels+1));out.tendency=h;
 for(std::size_t c=0;c<d.cells;++c){auto begin=c*d.levels;std::vector<Real> hm(h.air_mass.begin()+begin,h.air_mass.begin()+begin+d.levels),mass(d.air_mass_kg_m2.begin()+begin,d.air_mass_kg_m2.begin()+begin+d.levels);auto f=diagnose_vertical_mass_flux(hm,{b.begin(),b.end()},g);out.surface_pressure_pa_s[c]=f.surface_pressure_tendency_pa_s;out.maximum_continuity_residual_pa_s=std::max(out.maximum_continuity_residual_pa_s,std::abs(f.continuity_residual_pa_s));std::copy(f.interface_flux_kg_m2_s.begin(),f.interface_flux_kg_m2_s.end(),out.interface_mass_flux_kg_m2_s.begin()+static_cast<std::ptrdiff_t>(c*(d.levels+1)));
  auto apply=[&](const std::vector<Real>& state,const std::vector<Real>& horizontal){return vertical_scalar_rhs(state,mass,horizontal,f,scheme,limiter);}; std::vector<Real> x(d.levels),hx(d.levels);
  auto theta=apply({s.potential_temperature_mass_k_kg_m2.begin()+begin,s.potential_temperature_mass_k_kg_m2.begin()+begin+d.levels},{h.potential_temperature_mass.begin()+begin,h.potential_temperature_mass.begin()+begin+d.levels});auto tracer=apply({s.tracer_mass_kg_m2.begin()+begin,s.tracer_mass_kg_m2.begin()+begin+d.levels},{h.tracer_mass.begin()+begin,h.tracer_mass.begin()+begin+d.levels});for(std::size_t k=0;k<d.levels;++k){out.tendency.air_mass[begin+k]=f.target_air_mass_tendency_kg_m2_s[k];out.tendency.potential_temperature_mass[begin+k]=theta[k];out.tendency.tracer_mass[begin+k]=tracer[k];}
  for(int component=0;component<3;++component){for(std::size_t k=0;k<d.levels;++k){x[k]=s.horizontal_momentum_mass_kg_m_s[begin+k][component];hx[k]=h.momentum[begin+k][component];}auto rhs=apply(x,hx);for(std::size_t k=0;k<d.levels;++k)out.tendency.momentum[begin+k][component]=rhs[k];}
 } return out; }
}
