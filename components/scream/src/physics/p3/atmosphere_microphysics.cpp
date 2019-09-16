#include "share/scream_assert.hpp"
#include "physics/p3/scream_p3_interface.hpp"
#include "physics/p3/atmosphere_microphysics.hpp"
#include "share/simpleio/simple_io_mod.hpp"

namespace scream
{

// =========================================================================================
P3Microphysics::P3Microphysics (const Comm& comm,const ParameterList& /* params */)
  : m_p3_comm (comm)
{
/* Anything that can be initialized without grid information can be initialized here.
 * Like universal constants, table lookups, p3 options.
*/
}

// =========================================================================================
void P3Microphysics::set_grids(const std::shared_ptr<const GridsManager> grids_manager)
{

  using namespace units;
  // The units of mixing ratio Q are technically non-dimensional.
  // Nevertheless, for output reasons, we like to see 'kg/kg'.
  auto Q = kg/kg;
  auto Qdp = Q * Pa;
  Q.set_string("kg/kg");
  Qdp.set_string("kg/kg Pa");

  constexpr int NVL = 72;  /* TODO THIS NEEDS TO BE CHANGED TO A CONFIGURABLE */
  constexpr int QSZ =  9;  /* TODO THIS NEEDS TO BE CHANGED TO A CONFIGURABLE */

  auto grid = grids_manager->get_grid("Physics");
  const int num_dofs = grid->num_dofs();
  const int nc = 24; //num_dofs;

  auto VL = FieldTag::VerticalLevel;
  auto CO = FieldTag::Column;
  auto VR = FieldTag::Variable;
  auto TL = FieldTag::TimeLevel;

  FieldLayout scalar3d_layout { {CO,VL}, {nc,NVL} }; // Note that C++ and Fortran read array dimensions in reverse
  FieldLayout scalar3d_layout_int { {CO,VL}, {nc,NVL+1} }; // Note that C++ and Fortran read array dimensions in reverse
  FieldLayout vector3d_layout { {CO,VL, VR}, {nc, NVL, QSZ} };
  FieldLayout tracers_state_layout { {CO,TL,VR,VL}, {nc,2,4,NVL} };
  FieldLayout scalar_state_3d_mid_layout { {CO,TL,VL} , {nc,2,NVL}};
  FieldLayout q_forcing_layout  { {CO,VR,VL}, {nc,4,NVL} };

  // set requirements
  // Inputs from Dynamics
  m_required_fields.emplace("dp", scalar_state_3d_mid_layout,    Pa,             "Physics");
  m_required_fields.emplace("qdp", tracers_state_layout,        Qdp,          grid->name());  
  // Inputs from Physics 
  m_required_fields.emplace("ast",    scalar3d_layout,            A,           "Physics");
  m_required_fields.emplace("naai",   scalar3d_layout,   pow(kg,-1),           "Physics");
  m_required_fields.emplace("npccn",  scalar3d_layout,   pow(kg,-1)*pow(s,-1), "Physics");
  // Inputs that are STATE variables
  m_required_fields.emplace("pmid",  scalar3d_layout,    Pa,                   "Physics");
  m_required_fields.emplace("pdel",  scalar3d_layout,    Pa,                   "Physics");
  m_required_fields.emplace("zi",    scalar3d_layout_int, m,                   "Physics");
  m_required_fields.emplace("T",   scalar3d_layout,      K,         "Physics");
  m_required_fields.emplace("q",   vector3d_layout,      Q,         "Physics");
  m_required_fields.emplace("FQ", q_forcing_layout,      Q,      grid->name());

  // set computed
  // Inputs from Physics  NOTE: Making them "computed" fields now because we will initialize them in P3 for testing purposes.
  m_computed_fields.emplace("ast",    scalar3d_layout,            A,           "Physics");
  m_computed_fields.emplace("naai",   scalar3d_layout,   pow(kg,-1),           "Physics");
  m_computed_fields.emplace("npccn",  scalar3d_layout,   pow(kg,-1)*pow(s,-1), "Physics");
  m_computed_fields.emplace("pmid",  scalar3d_layout,    Pa,                   "Physics");
  m_computed_fields.emplace("pdel",  scalar3d_layout,    Pa,                   "Physics");
  m_computed_fields.emplace("zi",    scalar3d_layout_int, m,                   "Physics");
  // Fields changed by P3, not state variables 
  m_computed_fields.emplace("FQ", q_forcing_layout,      Q,      grid->name());
  // INPUT/OUTPUT from STATE
  m_computed_fields.emplace("T",   scalar3d_layout,      K,         "Physics");
  m_computed_fields.emplace("th",  scalar3d_layout,      K,         "Physics");
  m_computed_fields.emplace("q",   vector3d_layout,      Q,         "Physics");

}
// =========================================================================================
void P3Microphysics::initialize (const util::TimeStamp& t0)
{
  using namespace simpleio;
  m_current_ts = t0;
  //auto q_ptr = m_p3_fields_out.at("q").get_view().data();
  auto q_ptr = m_p3_fields_out.at("q").get_reshaped_view<Real[9][72][24]>();
  auto ast_ptr = m_p3_fields_out.at("ast").get_view().data();
  auto naai_ptr = m_p3_fields_out.at("naai").get_view().data();
  auto npccn_ptr = m_p3_fields_out.at("npccn").get_view().data();
  auto pmid_ptr = m_p3_fields_out.at("pmid").get_view().data();
  auto pdel_ptr = m_p3_fields_out.at("pdel").get_view().data();
  auto zi_ptr = m_p3_fields_out.at("zi").get_view().data();
  auto T_ptr = m_p3_fields_out.at("T").get_view().data();

  // Read input file information
  int rerr;
  int dlen[3];
  int flen  = 1728;
  int fleni = 1752;
  int time_lev = 1;
  char *filename;
  filename = "data/p3_convergence.nc";
  rerr = init_input(&filename);
  char *fieldname;
  char *units;
  // universal constants
  double uni_cnsts[12];
  dlen[0] = 12;
  fieldname = "universal_constants";
  rerr = readfield(&fieldname, time_lev, dlen, uni_cnsts); 
  // zi
  fieldname = "state%zi_in";
  dlen[0] = 24;
  dlen[1] = 73;
  rerr = readfield(&fieldname, time_lev, dlen, zi_ptr); 
  // T
  fieldname = "state%T_in";
  dlen[1] = 72;
  rerr = readfield(&fieldname, time_lev, dlen, T_ptr); 
  // PMID
  fieldname = "state%pmid_in";
  rerr = readfield(&fieldname, time_lev, dlen, pmid_ptr);
  // pdel
  fieldname = "state%pdel_in";
  rerr = readfield(&fieldname, time_lev, dlen, pdel_ptr); 
  // ast
  fieldname = "ast_in";
  rerr = readfield(&fieldname, time_lev, dlen, ast_ptr); 
  // naai 
  fieldname = "naai_in";
  rerr = readfield(&fieldname, time_lev, dlen, naai_ptr); 
  // npccn
  fieldname = "npccn_in";
  rerr = readfield(&fieldname, time_lev, dlen, npccn_ptr); 
  // Q is a bit more involved
  int ind;
  double *temp_ptr[9];
  char *q_list[9]= {"qv_in","cldliq_in","ice_in","numliq_in","numice_in","rain_in","numrain_in","qirim_in","rimvol_in"};
  for (int r=0;r<9;r++) {
    temp_ptr[r] = new double[flen];
    fieldname = q_list[r];
    rerr = readfield(&fieldname, time_lev, dlen, temp_ptr[r]);
    ind = 0;
    for (int k=0;k<72;k++) {
      for (int i=0;i<24;i++) {
        q_ptr(r,k,i) = temp_ptr[r][ind];
        ind++;
      } 
    }
  }

  rerr = finalize_input();

  std::string dimnames[] = {"column","level","ilevel"};
  const int ndims_output = sizeof(dimnames)/sizeof(dimnames[0]);
  int dimrng[] = {24,72,73};
  std::string vec_3d_state[] = {"column","level"};
  std::string vec_3d_inter[] = {"column","ilevel"};
  int werr, ftype;

  filename = "data/p3_output.nc";
  werr = init_output1(&filename,ndims_output,dimnames,dimrng);
  // Register fields for output
  fieldname = "th_out";
  units = "K";
  ftype = NC_REAL;
  werr = regfield(&fieldname,ftype,2,vec_3d_state,&units);
  werr = init_output2();

  p3_init_f90 (uni_cnsts);
}

// =========================================================================================
void P3Microphysics::run (const double dt)
{
  using namespace simpleio;
  auto q_ptr = m_p3_fields_out.at("q").get_view().data();
  auto FQ_ptr = m_p3_fields_out.at("FQ").get_view().data();
  auto dp_ptr = m_p3_fields_in.at("dp").get_view().data();
  auto qdp_ptr = m_p3_fields_in.at("qdp").get_view().data();
  auto ast_ptr = m_p3_fields_out.at("ast").get_view().data();
  auto naai_ptr = m_p3_fields_out.at("naai").get_view().data();
  auto npccn_ptr = m_p3_fields_out.at("npccn").get_view().data();
  auto pmid_ptr = m_p3_fields_out.at("pmid").get_view().data();
  auto pdel_ptr = m_p3_fields_out.at("pdel").get_view().data();
  auto zi_ptr = m_p3_fields_out.at("zi").get_view().data();
  auto T_ptr = m_p3_fields_out.at("T").get_view().data();
  auto th_ptr = m_p3_fields_out.at("th").get_view().data();

  p3_main_f90 (dt,q_ptr,FQ_ptr,qdp_ptr,T_ptr,zi_ptr,pmid_ptr,pdel_ptr,
     ast_ptr, naai_ptr, npccn_ptr, th_ptr);

  m_current_ts += dt;
  m_p3_fields_out.at("q").get_header().get_tracking().update_time_stamp(m_current_ts);
  m_p3_fields_out.at("FQ").get_header().get_tracking().update_time_stamp(m_current_ts);
  m_p3_fields_out.at("T").get_header().get_tracking().update_time_stamp(m_current_ts);
  m_p3_fields_out.at("th").get_header().get_tracking().update_time_stamp(m_current_ts);

  char *fieldname;
  int ftype, werr, ind;
  int dim3d[]  = {24,72};
  int dim3di[] = {24,73};
//  double tmp2d[72][24];

  fieldname = "th_out";
//  ind = 0;
//  for (int j=0;j<72;j++) {
//    for (int i=0;i<24;i++) {
//      tmp2d[j][i] = th_ptr[ind];
//      ind++;
//    }
//  }
  werr = writefield(&fieldname,1,dim3d,th_ptr);

}
// =========================================================================================
void P3Microphysics::finalize()
{
  using namespace simpleio;
  p3_finalize_f90 ();
  int werr = finalize_output();
}
// =========================================================================================

void P3Microphysics::register_fields (FieldRepository<Real, device_type>& field_repo) const {
  for (auto& fid : m_required_fields) {
    field_repo.register_field(fid);
  }
  for (auto& fid : m_computed_fields) {
    field_repo.register_field(fid);
  }
}

void P3Microphysics::set_required_field_impl (const Field<const Real, device_type>& f) {
  // Store a copy of the field. We need this in order to do some tracking checks
  // at the beginning of the run call. Other than that, there would be really
  // no need to store a scream field here; we could simply set the view ptr
  // in the Homme's view, and be done with it.
  m_p3_fields_in.emplace(f.get_header().get_identifier().name(),f);

  // Add myself as customer to the field
  f.get_header_ptr()->get_tracking().add_customer(weak_from_this());
}

void P3Microphysics::set_computed_field_impl (const Field<      Real, device_type>& f) {
  // Store a copy of the field. We need this in order to do some tracking updates
  // at the end of the run call. Other than that, there would be really
  // no need to store a scream field here; we could simply set the view ptr
  // in the Homme's view, and be done with it.
  m_p3_fields_out.emplace(f.get_header().get_identifier().name(),f);

  // Add myself as provider for the field
  f.get_header_ptr()->get_tracking().add_provider(weak_from_this());
}

} // namespace scream
