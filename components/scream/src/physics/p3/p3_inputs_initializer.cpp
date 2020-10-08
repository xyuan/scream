#include "physics/p3/p3_inputs_initializer.hpp"
#include "physics/p3/scream_p3_interface.hpp"

#include <array>

namespace scream
{

const int INPUT_SIZE = 44;

//Layout options are set as an int 
//to be passed into the GridOpts struct
const int SCALAR_3D_MID = 0;
const int SCALAR_3D_INT = 1;
const int VECTOR_3D_MID = 2;
const int TRACERS = 3;
const int LINEAR = 4;

void P3InputsInitializer::add_field (const field_type &f)
{
p3_inputs = {"q","T","FQ","ast","ni_activated",
		"nc_nuceat_tend","pmid","dp","zi", "qc", "nc", 
		"qr" , "nr", "qi", "qm", "ni", "bm", "qv", "th",
		"inv_qc_relvar", "cld_frac_i", "cld_frac_l", "cld_frac_r", 
		"pres", "dz", "dpres", "exner", "liq_ice_exchange", 
		"vap_liq_exchange", "vap_ice_exchange", "mu_c", "lamc",
		"cmeiout", "precip_liq_surf", "precip_ice_surf", "diag_effc",
		"diag_effi", "rho_qi", "precip_total_tend", "nevapr",
		"qr_evap_tend", "precip_liq_flux", "precip_ice_flux",
		"col_location"};

  const auto& id = f.get_header().get_identifier();
  
  m_fields.emplace(id.name(),f);
  m_fields_id.insert(id);
}

// =========================================================================================
void P3InputsInitializer::initialize_fields ()
{
  // Safety check: if we're asked to init anything at all,
  // then we should have been asked to init 7 fields.
//  int count = 0;
//  count += m_fields.count("q");
//  count += m_fields.count("T");
//  count += m_fields.count("ast");
//  count += m_fields.count("ni_activated");
//  count += m_fields.count("nc_nuceat_tend");
//  count += m_fields.count("pmid");
//  count += m_fields.count("dp");
//  count += m_fields.count("zi");
//  
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

  // Get device views
  auto d_q     = m_fields.at("q").get_view();
  auto d_T     = m_fields.at("T").get_view();
  auto d_ast   = m_fields.at("ast").get_view();
  auto d_ni_activated  = m_fields.at("ni_activated").get_view();
  auto d_nc_nuceat_tend = m_fields.at("nc_nuceat_tend").get_view();
  auto d_pmid  = m_fields.at("pmid").get_view();
  auto d_dpres  = m_fields.at("dp").get_view();
  auto d_zi    = m_fields.at("zi").get_view();

  // Create host mirrors
  auto h_q     = Kokkos::create_mirror_view(d_q);
  auto h_T     = Kokkos::create_mirror_view(d_T);
  auto h_ast   = Kokkos::create_mirror_view(d_ast);
  auto h_ni_activated  = Kokkos::create_mirror_view(d_ni_activated);
  auto h_nc_nuceat_tend = Kokkos::create_mirror_view(d_nc_nuceat_tend);
  auto h_pmid  = Kokkos::create_mirror_view(d_pmid);
  auto h_dpres  = Kokkos::create_mirror_view(d_dpres);
  auto h_zi    = Kokkos::create_mirror_view(d_zi);

  // Get host mirros' raw pointers
  auto q     = h_q.data();
  auto T     = h_T.data();
  auto ast   = h_ast.data();
  auto ni_activated  = h_ni_activated.data();
  auto nc_nuceat_tend = h_nc_nuceat_tend.data();
  auto pmid  = h_pmid.data();
  auto dpres  = h_dpres.data();
  auto zi    = h_zi.data();

  // Call f90 routine
  p3_standalone_init_f90 (q, T, zi, pmid, dpres, ast, ni_activated, nc_nuceat_tend);

  // Deep copy back to device
  Kokkos::deep_copy(d_q,h_q);
  Kokkos::deep_copy(d_T,h_T);
  Kokkos::deep_copy(d_ast,h_ast);
  Kokkos::deep_copy(d_ni_activated,h_ni_activated);
  Kokkos::deep_copy(d_nc_nuceat_tend,h_nc_nuceat_tend);
  Kokkos::deep_copy(d_pmid,h_pmid);
  Kokkos::deep_copy(d_dpres,h_dpres);
  Kokkos::deep_copy(d_zi,h_zi);

  // If we are in charge of init-ing FQ as well, init it to 0.
  if (m_fields.count("FQ")==1) {
    // Init FQ to 0
    auto d_FQ = m_fields.at("FQ").get_view();
    Kokkos::deep_copy(d_FQ,Real(0));
  }
}

} // namespace scream
