#include "physics/p3/p3_inputs_initializer.hpp"
#include "physics/p3/scream_p3_interface.hpp"

#include <array>
#include <unordered_map>

namespace scream
{

using namespace std;

vector<string> p3_inputs = {"q","T","FQ","ast","ni_activated",
		"nc_nuceat_tend","pmid","dp","zi", "qc", "nc", 
		"qr" , "nr", "qi", "qm", "ni", "bm", "qv", "th",
		"inv_qc_relvar", "cld_frac_i", "cld_frac_l", "cld_frac_r", 
		"pres", "dz", "dpres", "exner", "liq_ice_exchange", 
		"vap_liq_exchange", "vap_ice_exchange", "mu_c", "lamc",
		"cmeiout", "precip_liq_surf", "precip_ice_surf", "diag_effc",
		"diag_effi", "rho_qi", "precip_total_tend", "nevapr",
		"qr_evap_tend", "precip_liq_flux", "precip_ice_flux",
		"col_location"};


const int INPUT_SIZE = 44;

//Layout options are set as an int 
//to be passed into the GridOpts struct
const int SCALAR_3D_MID = 0;
const int SCALAR_3D_INT = 1;
const int VECTOR_3D_MID = 2;
const int TRACERS = 3;
const int LINEAR = 4;

using namespace scream;
using namespace ekat;
using namespace units;

auto Q = kg/kg;
auto nondim = m/m;




struct GridOpts{
  string name;
  bool isOut;
  const Units* unit;
  int field_idx;
};


unordered_map<string, GridOpts> opt_map;

//Initializes struct GridOpts fields
void set_grid_opts_helper(GridOpts O, string n, bool out, const Units* unit, int field_idx
                          ){
  
  O.name = n;
  O.isOut = out;
  O.field_idx = field_idx;
  O.unit = unit;
  opt_map.insert({O.name, O});
}

void set_grid_opts(){

  GridOpts q; 
  GridOpts T;
  GridOpts FQ;
  GridOpts ast;
  GridOpts ni_activated;
  GridOpts nc_nuceat_tend;
  GridOpts pmid;
  GridOpts dp;
  GridOpts zi;
  GridOpts qc;
  GridOpts nc;
  GridOpts qr;
  GridOpts nr;
  GridOpts qi;
  GridOpts qm;
  GridOpts ni;
  GridOpts bm;
  GridOpts qv;
  GridOpts th;
  GridOpts inv_qc_relvar;
  GridOpts cld_frac_i;
  GridOpts cld_frac_l;
  GridOpts cld_frac_r;
  GridOpts pres;
  GridOpts dz;
  GridOpts dpres;
  GridOpts exner;
  GridOpts liq_ice_exchange;
  GridOpts vap_liq_exchange;
  GridOpts vap_ice_exchange;
  GridOpts mu_c;
  GridOpts lamc;
  GridOpts cmeiout;
  GridOpts precip_liq_surf;
  GridOpts precip_ice_surf;
  GridOpts diag_effc;
  GridOpts diag_effi;
  GridOpts rho_qi;
  GridOpts precip_total_tend;
  GridOpts nevapr;
  GridOpts qr_evap_tend;
  GridOpts precip_liq_flux;
  GridOpts precip_ice_flux;
  GridOpts col_location;
  
  
  set_grid_opts_helper(ast, "ast", true, &Q, SCALAR_3D_MID);
  set_grid_opts_helper(ni_activated, "ni_activated", true, &(1/kg), SCALAR_3D_MID);
  set_grid_opts_helper(nc_nuceat_tend, "nc_nuceat_tend", true, &(1/(kg*s)), SCALAR_3D_MID);
  set_grid_opts_helper(pmid, "pmid", true, &Pa, SCALAR_3D_MID);
  set_grid_opts_helper(dp, "dp", true, &Pa, SCALAR_3D_MID);
  set_grid_opts_helper(zi, "zi", true, &Q, SCALAR_3D_MID);
  set_grid_opts_helper(qc, "qc", true, &Pa, TRACERS);
  set_grid_opts_helper(nc, "nc", true, &Pa, TRACERS);
  set_grid_opts_helper(qr, "qr", true, &Pa, TRACERS);
  set_grid_opts_helper(nr, "nr", true, &Pa, TRACERS);
  set_grid_opts_helper(qi, "qi", true, &Pa, TRACERS);
  set_grid_opts_helper(qm, "qm", true, &Pa, TRACERS);
  set_grid_opts_helper(ni, "ni", true, &Pa, TRACERS);
  set_grid_opts_helper(bm, "bm", true, &Pa, TRACERS);
  set_grid_opts_helper(qv, "qv", true, &Pa, TRACERS);
  set_grid_opts_helper(th, "th", true, &Pa, TRACERS);
  set_grid_opts_helper(inv_qc_relvar, "inv_qc_relvar", true, &Pa, SCALAR_3D_INT);
  set_grid_opts_helper(cld_frac_i, "cld_frac_i", true, &Pa, SCALAR_3D_INT);
  set_grid_opts_helper(cld_frac_l, "cld_frac_l", true, &Pa, SCALAR_3D_INT);
  set_grid_opts_helper(cld_frac_r, "cld_frac_r", true, &Pa, SCALAR_3D_INT);
  set_grid_opts_helper(pres, "pres", true, &Pa, SCALAR_3D_INT);
  set_grid_opts_helper(dz, "dz", true, &Pa, SCALAR_3D_INT);
  set_grid_opts_helper(dpres, "dpres", true, &Pa, SCALAR_3D_INT);
  set_grid_opts_helper(exner, "exner", true, &Pa, SCALAR_3D_INT);
  set_grid_opts_helper(liq_ice_exchange, "liq_ice_exchange", true, &Pa, SCALAR_3D_INT);
  set_grid_opts_helper(vap_liq_exchange, "vap_liq_exchange", true, &Pa, SCALAR_3D_INT);
  set_grid_opts_helper(vap_ice_exchange, "vap_ice_exchange", true, &Pa, SCALAR_3D_INT);
  set_grid_opts_helper(col_location, "col_location", true, &Pa, SCALAR_3D_INT);
  set_grid_opts_helper(mu_c, "mu_c", true, &Pa, SCALAR_3D_INT);
  set_grid_opts_helper(lamc, "lamc", true, &Pa, SCALAR_3D_INT);
  set_grid_opts_helper(cmeiout, "cmeiout", true, &Pa, SCALAR_3D_INT);
  set_grid_opts_helper(precip_liq_surf, "precip_liq_surf", true, &Pa, SCALAR_3D_INT);
  set_grid_opts_helper(precip_ice_surf, "precip_ice_surf", true, &Pa, SCALAR_3D_INT);
  set_grid_opts_helper(diag_effc, "diag_effc", true, &Pa, SCALAR_3D_INT);
  set_grid_opts_helper(diag_effi, "diag_effi", true, &Pa, SCALAR_3D_INT);
  set_grid_opts_helper(rho_qi, "rho_qi", true, &Pa, SCALAR_3D_INT);
  set_grid_opts_helper(precip_total_tend, "precip_total_tend", true, &Pa, SCALAR_3D_INT);
  set_grid_opts_helper(nevapr, "nevapr", true, &Pa, SCALAR_3D_INT);
  set_grid_opts_helper(qr_evap_tend, "qr_evap_tend", true, &Pa, SCALAR_3D_INT);
  set_grid_opts_helper(precip_liq_flux, "precip_liq_flux", true, &Pa, SCALAR_3D_INT);
  set_grid_opts_helper(precip_ice_flux, "precip_ice_flux", true, &Pa, SCALAR_3D_INT);
  set_grid_opts_helper(q, "q", true, &Q, TRACERS);
  set_grid_opts_helper(T, "t", true, &K, SCALAR_3D_MID);
  set_grid_opts_helper(FQ, "FQ", true, &Q, VECTOR_3D_MID);


  
  
  }
  void P3InputsInitializer::add_field (const field_type &f)
  {
    const auto& id = f.get_header().get_identifier();
    
    m_fields.emplace(id.name(),f);
    m_fields_id.insert(id);
  }
  
  // =========================================================================================
  void P3InputsInitializer::initialize_fields ()
{
  // Safety check: if we're asked to init anything at all,
  // then we should have been asked to init 7 fields.
//  check if p3 inputs have been registered as fields
  int count = 0;
  for (int j = 0; j < p3_inputs.size(); j++){
    count += m_fields.count(p3_inputs[j]);
  } 
  if (count==0) {
    return;
  }

  EKAT_REQUIRE_MSG (count==INPUT_SIZE,
    "Error! P3InputsInitializer is expected to init 'q','T','ast','ni_activated','nc_nuceat_tend','pmid','dp','zi'.\n"
    "       Only " + std::to_string(count) + " of those have been found.\n"
    "       Please, check the atmosphere processes you are using,"
    "       and make sure they agree on who's initializing each field.\n");

  //// Get device views
  //auto d_q     = m_fields.at("q").get_view();
  //auto d_T     = m_fields.at("T").get_view();
  //auto d_ast   = m_fields.at("ast").get_view();
  //auto d_ni_activated  = m_fields.at("ni_activated").get_view();
  //auto d_nc_nuceat_tend = m_fields.at("nc_nuceat_tend").get_view();
  //auto d_pmid  = m_fields.at("pmid").get_view();
  //auto d_dpres  = m_fields.at("dp").get_view();
  //auto d_zi    = m_fields.at("zi").get_view();

  //// Create host mirrors
  //auto h_q     = Kokkos::create_mirror_view(d_q);
  //auto h_T     = Kokkos::create_mirror_view(d_T);
  //auto h_ast   = Kokkos::create_mirror_view(d_ast);
  //auto h_ni_activated  = Kokkos::create_mirror_view(d_ni_activated);
  //auto h_nc_nuceat_tend = Kokkos::create_mirror_view(d_nc_nuceat_tend);
  //auto h_pmid  = Kokkos::create_mirror_view(d_pmid);
  //auto h_dpres  = Kokkos::create_mirror_view(d_dpres);
  //auto h_zi    = Kokkos::create_mirror_view(d_zi);

  //// Get host mirros' raw pointers
  //auto q     = h_q.data();
  //auto T     = h_T.data();
  //auto ast   = h_ast.data();
  //auto ni_activated  = h_ni_activated.data();
  //auto nc_nuceat_tend = h_nc_nuceat_tend.data();
  //auto pmid  = h_pmid.data();
  //auto dpres  = h_dpres.data();
  //auto zi    = h_zi.data();

  // Call f90 routine
//  p3_standalone_init_f90 (q, T, zi, pmid, dpres, ast, ni_activated, nc_nuceat_tend);

  // Deep copy back to device
  //Kokkos::deep_copy(d_q,h_q);
  //Kokkos::deep_copy(d_T,h_T);
  //Kokkos::deep_copy(d_ast,h_ast);
  //Kokkos::deep_copy(d_ni_activated,h_ni_activated);
  //Kokkos::deep_copy(d_nc_nuceat_tend,h_nc_nuceat_tend);
  //Kokkos::deep_copy(d_pmid,h_pmid);
  //Kokkos::deep_copy(d_dpres,h_dpres);
  //Kokkos::deep_copy(d_zi,h_zi);


  for (int i = 0; i < p3_inputs.size(); i++){
    //Get and store device view using input name
    Kokkos::View<Real*, Kokkos::LayoutRight, HostDevice> d_v = m_fields.at(p3_inputs[i]).get_view();
    //Create and store host mirrors using device views
    Kokkos::View<scream::Real*, Kokkos::LayoutRight, HostDevice> h_m = Kokkos::create_mirror_view(d_v);
    //Create and store host mirrors raw pointers
    Real* r_p = h_m.data();
    //Deep copy back to device
    Kokkos::deep_copy(d_v, h_m);
  }

  // If we are in charge of init-ing FQ as well, init it to 0.
  if (m_fields.count("FQ")==1) {
    // Init FQ to 0
    auto d_FQ = m_fields.at("FQ").get_view();
    Kokkos::deep_copy(d_FQ,Real(0));
  }
}

} // namespace scream
