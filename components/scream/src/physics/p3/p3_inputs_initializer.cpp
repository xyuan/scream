#include "physics/p3/p3_inputs_initializer.hpp"
#include "physics/share/physics_constants.hpp"
#include "physics/p3/p3_main_impl.hpp"
#include "share/io/scorpio_input.hpp"
#include "share/io/scream_scorpio_interface.hpp"

#include "ekat/util/ekat_file_utils.hpp"
#include <array>
#include <fstream>

namespace scream
{

void P3InputsInitializer::add_field (const field_type &f)
{
  const auto& id = f.get_header().get_identifier();
  
  m_fields.emplace(id.name(),f);
  m_fields_id.insert(id);
}

void P3InputsInitializer::
add_field (const field_type &f, const field_type& f_ref,
           const remapper_ptr_type& remapper)
{
  if (m_remapper) {
    // Sanity check
    EKAT_REQUIRE_MSG (m_remapper->get_src_grid()->name()==remapper->get_src_grid()->name(),
      "Error! A remapper was already set in P3InputsInitializer, but its src grid differs from"
      "       the grid of the input remapper of this call.\n");
  } else {
    m_remapper = remapper;
    m_remapper->registration_begins();
  }

  const auto& id = f.get_header().get_identifier();
  const auto& id_ref = f_ref.get_header().get_identifier();

  // To the AD, we only expose the fact that we init f_ref...
  m_fields_id.insert(id_ref);

  // ...but P3 only knows how to init f...
  m_fields.emplace(id.name(),f);

  // ...hence, we remap to f_ref.
  m_remapper->register_field(f, f_ref);
}


// =========================================================================================
void P3InputsInitializer::initialize_fields ()
{
  using namespace p3;
  using input_type = AtmosphereInput;
  printf("--------------------------- Start p3 initializer ----------------------\n");
  // To simplify the initializer we first define all the fields we expect to have to initialize.
  std::vector<std::string> fields_to_init;
  fields_to_init.push_back("T_atm");
  fields_to_init.push_back("ast");
  fields_to_init.push_back("ni_activated");
  fields_to_init.push_back("nc_nuceat_tend");
  fields_to_init.push_back("pmid");
  fields_to_init.push_back("dp");
  fields_to_init.push_back("zi");
  fields_to_init.push_back("qv_prev");
  fields_to_init.push_back("T_prev");
  fields_to_init.push_back("qv");
  fields_to_init.push_back("qc");
  fields_to_init.push_back("qr");
  fields_to_init.push_back("qi");
  fields_to_init.push_back("qm");
  fields_to_init.push_back("nc");
  fields_to_init.push_back("nr");
  fields_to_init.push_back("ni");
  fields_to_init.push_back("bm");
  fields_to_init.push_back("nccn_prescribed");
  fields_to_init.push_back("inv_qc_relvar");
  // TODO: Delete eventually, should instead be set in run interface: 
  fields_to_init.push_back("th_atm");  
  fields_to_init.push_back("dz");  
  fields_to_init.push_back("exner");  
  fields_to_init.push_back("cld_frac_l");  
  fields_to_init.push_back("cld_frac_i");  
  fields_to_init.push_back("cld_frac_r");  
  // Safety check: if we're asked to init anything at all,
  // then we should have been asked to init 20 fields.
  int count = 0;
  std::string list_of_fields = "";
  for (auto name : fields_to_init)
  {
    list_of_fields += name;
    list_of_fields += ", ";
    count += m_fields.count(name);
  }
  
  if (count==0) {
    return;
  }

  EKAT_REQUIRE_MSG (count==fields_to_init.size(),
//    "Error! P3InputsInitializer is expected to init 'q','T','ast','ni_activated','nc_nuceat_tend','pmid','dp','zi','qv_prev','T_prev'.\n"
    "Error! P3InputsInitializer is expected to init " + std::to_string(fields_to_init.size()) + " fields:\n"
    "       " + list_of_fields + "\n"
    "       Instead found " + std::to_string(count) + " fields.\n"
    "       Please, check the atmosphere processes you are using,\n"
    "       and make sure they agree on who's initializing each field.\n");

  // Get views
  using P3F             = Functions<Real, DefaultDevice>;
  using Spack           = typename P3F::Spack;
  using value_type      = Real;
  using device_type     = DefaultDevice;
  using Pack            = ekat::Pack<Real,Spack::n>;
  using DeviceView      = ekat::Unmanaged<typename KokkosTypes<DefaultDevice>::template view<Pack**> >;
  using HostView        = ekat::Unmanaged<typename KokkosTypes<HostDevice>::template view<Pack**> >;

  std::map<std::string,DeviceView> device_views;
  std::map<std::string,HostView> host_mirrors;
  for (auto name : fields_to_init)
  {
    // Get device views
    //DeviceView d_view =  m_fields.at(name).get_reshaped_view<Real**>();
    auto d_view_tmp = m_fields.at(name).get_reshaped_view<Pack**>();
    DeviceView d_view =  m_fields.at(name).get_reshaped_view<Pack**>();
    EKAT_REQUIRE_MSG(typeid(d_view_tmp).name()==typeid(d_view).name(),"device types are not the same,");// + typeid(d_view_tmp).name().c_str() + " vs. " + typeid(d_view).name().c_str());
    device_views.emplace(name,d_view);
    // Create host mirrors
    auto h_view     = Kokkos::create_mirror_view(d_view);
    auto h_view_tmp = Kokkos::create_mirror_view(d_view_tmp);
    EKAT_REQUIRE_MSG(typeid(h_view_tmp).name()==typeid(h_view).name(),"host types are not the same,");// + typeid(h_view_tmp).name().c_str() + " vs. " + typeid(h_view).name().c_str());
    host_mirrors.emplace(name,h_view);//Kokkos::create_mirror_view(d_view));
  }

  // Create device views
  auto d_T_atm           = m_fields.at("T_atm").get_reshaped_view<Pack**>();
  auto d_ast             = m_fields.at("ast").get_reshaped_view<Pack**>();
  auto d_ni_activated    = m_fields.at("ni_activated").get_reshaped_view<Pack**>();
  auto d_nc_nuceat_tend  = m_fields.at("nc_nuceat_tend").get_reshaped_view<Pack**>();
  auto d_pmid            = m_fields.at("pmid").get_reshaped_view<Pack**>();
  auto d_dp              = m_fields.at("dp").get_reshaped_view<Pack**>();
  auto d_zi              = m_fields.at("zi").get_reshaped_view<Pack**>();
  auto d_qv_prev         = m_fields.at("qv_prev").get_reshaped_view<Pack**>();
  auto d_T_prev          = m_fields.at("T_prev").get_reshaped_view<Pack**>();
  auto d_qv              = m_fields.at("qv").get_reshaped_view<Pack**>();
  auto d_qc              = m_fields.at("qc").get_reshaped_view<Pack**>();
  auto d_qr              = m_fields.at("qr").get_reshaped_view<Pack**>();
  auto d_qi              = m_fields.at("qi").get_reshaped_view<Pack**>();
  auto d_qm              = m_fields.at("qm").get_reshaped_view<Pack**>();
  auto d_nc              = m_fields.at("nc").get_reshaped_view<Pack**>();
  auto d_nr              = m_fields.at("nr").get_reshaped_view<Pack**>();
  auto d_ni              = m_fields.at("ni").get_reshaped_view<Pack**>();
  auto d_bm              = m_fields.at("bm").get_reshaped_view<Pack**>();
  auto d_nccn_prescribed = m_fields.at("nccn_prescribed").get_reshaped_view<Pack**>();
  auto d_inv_qc_relvar   = m_fields.at("inv_qc_relvar").get_reshaped_view<Pack**>();
  // TODO: Delete eventually, should instead be set in run interface: 
  auto d_th_atm          = m_fields.at("th_atm").get_reshaped_view<Pack**>();  
  auto d_dz              = m_fields.at("dz").get_reshaped_view<Pack**>();  
  auto d_exner           = m_fields.at("exner").get_reshaped_view<Pack**>();  
  auto d_cld_frac_l      = m_fields.at("cld_frac_l").get_reshaped_view<Pack**>();  
  auto d_cld_frac_i      = m_fields.at("cld_frac_i").get_reshaped_view<Pack**>();  
  auto d_cld_frac_r      = m_fields.at("cld_frac_r").get_reshaped_view<Pack**>();
  // Create Host Views
  auto h_T_atm            = Kokkos::create_mirror_view(d_T_atm          ); 
  auto h_ast              = Kokkos::create_mirror_view(d_ast            ); 
  auto h_ni_activated     = Kokkos::create_mirror_view(d_ni_activated   ); 
  auto h_nc_nuceat_tend   = Kokkos::create_mirror_view(d_nc_nuceat_tend ); 
  auto h_pmid             = Kokkos::create_mirror_view(d_pmid           ); 
  auto h_dp               = Kokkos::create_mirror_view(d_dp             ); 
  auto h_zi               = Kokkos::create_mirror_view(d_zi             ); 
  auto h_qv_prev          = Kokkos::create_mirror_view(d_qv_prev        ); 
  auto h_T_prev           = Kokkos::create_mirror_view(d_T_prev         ); 
  auto h_qv               = Kokkos::create_mirror_view(d_qv             ); 
  auto h_qc               = Kokkos::create_mirror_view(d_qc             ); 
  auto h_qr               = Kokkos::create_mirror_view(d_qr             ); 
  auto h_qi               = Kokkos::create_mirror_view(d_qi             ); 
  auto h_qm               = Kokkos::create_mirror_view(d_qm             ); 
  auto h_nc               = Kokkos::create_mirror_view(d_nc             ); 
  auto h_nr               = Kokkos::create_mirror_view(d_nr             ); 
  auto h_ni               = Kokkos::create_mirror_view(d_ni             ); 
  auto h_bm               = Kokkos::create_mirror_view(d_bm             ); 
  auto h_nccn_prescribed  = Kokkos::create_mirror_view(d_nccn_prescribed); 
  auto h_inv_qc_relvar    = Kokkos::create_mirror_view(d_inv_qc_relvar  ); 
  // TODO: Delete eventually, should instead be set in run interface: 
  auto h_th_atm           = Kokkos::create_mirror_view(d_th_atm         ); 
  auto h_dz               = Kokkos::create_mirror_view(d_dz             ); 
  auto h_exner            = Kokkos::create_mirror_view(d_exner          ); 
  auto h_cld_frac_l       = Kokkos::create_mirror_view(d_cld_frac_l     ); 
  auto h_cld_frac_i       = Kokkos::create_mirror_view(d_cld_frac_i     ); 
  auto h_cld_frac_r       = Kokkos::create_mirror_view(d_cld_frac_r     ); 
  // Initalize from text file 
  std::ifstream fid("p3_init_vals.txt", std::ifstream::in);
  std::string tmp_line;
  int icol_in_max = 0;
  while(getline(fid,tmp_line))
  {
    std::stringstream s(tmp_line);
    std::string field;
    std::vector<Real> field_vals;
    while (getline(s,field,' '))
    {
      field_vals.push_back(std::stod(field));
    }
    int icol  = (int)field_vals[0];
    int ipack = (int)field_vals[1] / Spack::n;
    int ivec  = (int)field_vals[1] % Spack::n;
    icol_in_max = std::max(icol,icol_in_max);
    int cnt = 1;
    cnt++;
    h_qv(icol,ipack)[ivec]=field_vals[cnt];
    cnt++;
    h_th_atm(icol,ipack)[ivec]=field_vals[cnt];
    cnt++;
    h_pmid(icol,ipack)[ivec]=field_vals[cnt];
    cnt++;
    h_dz(icol,ipack)[ivec]=field_vals[cnt];
    cnt++;
    h_nc_nuceat_tend(icol,ipack)[ivec]=field_vals[cnt];
    cnt++;
    h_nccn_prescribed(icol,ipack)[ivec]=field_vals[cnt];
    cnt++;
    h_ni_activated(icol,ipack)[ivec]=field_vals[cnt];
    cnt++;
    h_inv_qc_relvar(icol,ipack)[ivec]=field_vals[cnt];
    cnt++;
    h_qc(icol,ipack)[ivec]=field_vals[cnt];
    cnt++;
    h_nc(icol,ipack)[ivec]=field_vals[cnt];
    cnt++;
    h_qr(icol,ipack)[ivec]=field_vals[cnt];
    cnt++;
    h_nr(icol,ipack)[ivec]=field_vals[cnt];
    cnt++;
    h_qi(icol,ipack)[ivec]=field_vals[cnt];
    cnt++;
    h_ni(icol,ipack)[ivec]=field_vals[cnt];
    cnt++;
    h_qm(icol,ipack)[ivec]=field_vals[cnt];
    cnt++;
    h_bm(icol,ipack)[ivec]=field_vals[cnt];
    cnt++;
    //h_precip_liq_surf(icol,ipack)[ivec]=field_vals[cnt];    //
    cnt++;
    //h_precip_ice_surf(icol,ipack)[ivec]=field_vals[cnt];    //
    cnt++;
    //h_diag_eff_radius_qc(icol,ipack)[ivec]=field_vals[cnt]; //
    cnt++;
    //h_diag_eff_radius_qi(icol,ipack)[ivec]=field_vals[cnt]; //
    cnt++;
    //h_rho_qi(icol,ipack)[ivec]=field_vals[cnt]; //
    cnt++;
    h_dp(icol,ipack)[ivec]=field_vals[cnt];
    cnt++;
    h_exner(icol,ipack)[ivec]=field_vals[cnt];
    cnt++;
    //h_qv2qi_depos_tend(icol,ipack)[ivec]=field_vals[cnt]; //
    cnt++;
    //h_precip_total_tend(icol,ipack)[ivec]=field_vals[cnt]; //
    cnt++;
    //h_nevapr(icol,ipack)[ivec]=field_vals[cnt]; //
    cnt++;
    //h_qr_evap_tend(icol,ipack)[ivec]=field_vals[cnt]; //
    cnt++;
    //h_precip_liq_flux(icol,ipack)[ivec]=field_vals[cnt]; //
    cnt++;
    //h_precip_ice_flux(icol,ipack)[ivec]=field_vals[cnt]; //
    cnt++;
    h_cld_frac_r(icol,ipack)[ivec]=field_vals[cnt];
    cnt++;
    h_cld_frac_l(icol,ipack)[ivec]=field_vals[cnt];
    cnt++;
    h_cld_frac_i(icol,ipack)[ivec]=field_vals[cnt];
    cnt++;
    //h_mu_c(icol,ipack)[ivec]=field_vals[cnt]; //
    cnt++;
    //h_lamc(icol,ipack)[ivec]=field_vals[cnt]; //
    cnt++;
    //h_liq_ice_exchange(icol,ipack)[ivec]=field_vals[cnt]; //
    cnt++;
    //h_vap_liq_exchange(icol,ipack)[ivec]=field_vals[cnt]; //
    cnt++;
    //h_vap_ice_exchange(icol,ipack)[ivec]=field_vals[cnt]; //
    cnt++;
    h_qv_prev(icol,ipack)[ivec]=field_vals[cnt]; 
    cnt++;
    h_T_prev(icol,ipack)[ivec]=field_vals[cnt];
  }
  // For now use dummy values copied from `p3_ic_cases.cpp`
  using consts          = scream::physics::Constants<Real>;
  auto temp = m_fields.at("qc");
  auto mdims = temp.get_header().get_identifier().get_layout();
  Int ncol = mdims.dim(0); 
  Int nk   = mdims.dim(1);

  for (int icol_i = icol_in_max;icol_i<ncol;icol_i++)
  {
    for (int k = 0;k<nk;k++)
    {
    int icol  = icol_i % icol_in_max;
    int ipack = k / Spack::n;
    int ivec  = k % Spack::n;
    h_qv(icol_i,ipack)[ivec]              = h_qv(icol,ipack)[ivec]             ;
    h_th_atm(icol_i,ipack)[ivec]          = h_th_atm(icol,ipack)[ivec]         ;
    h_pmid(icol_i,ipack)[ivec]            = h_pmid(icol,ipack)[ivec]           ;
    h_dz(icol_i,ipack)[ivec]              = h_dz(icol,ipack)[ivec]             ;
    h_nc_nuceat_tend(icol_i,ipack)[ivec]  = h_nc_nuceat_tend(icol,ipack)[ivec] ;
    h_nccn_prescribed(icol_i,ipack)[ivec] = h_nccn_prescribed(icol,ipack)[ivec];
    h_ni_activated(icol_i,ipack)[ivec]    = h_ni_activated(icol,ipack)[ivec]   ;
    h_inv_qc_relvar(icol_i,ipack)[ivec]   = h_inv_qc_relvar(icol,ipack)[ivec]  ;
    h_qc(icol_i,ipack)[ivec]              = h_qc(icol,ipack)[ivec]             ;
    h_nc(icol_i,ipack)[ivec]              = h_nc(icol,ipack)[ivec]             ;
    h_qr(icol_i,ipack)[ivec]              = h_qr(icol,ipack)[ivec]             ;
    h_nr(icol_i,ipack)[ivec]              = h_nr(icol,ipack)[ivec]             ;
    h_qi(icol_i,ipack)[ivec]              = h_qi(icol,ipack)[ivec]             ;
    h_ni(icol_i,ipack)[ivec]              = h_ni(icol,ipack)[ivec]             ;
    h_qm(icol_i,ipack)[ivec]              = h_qm(icol,ipack)[ivec]             ;
    h_bm(icol_i,ipack)[ivec]              = h_bm(icol,ipack)[ivec]             ;
    h_dp(icol_i,ipack)[ivec]              = h_dp(icol,ipack)[ivec]             ;
    h_exner(icol_i,ipack)[ivec]           = h_exner(icol,ipack)[ivec]          ;
    h_cld_frac_r(icol_i,ipack)[ivec]      = h_cld_frac_r(icol,ipack)[ivec]     ;
    h_cld_frac_l(icol_i,ipack)[ivec]      = h_cld_frac_l(icol,ipack)[ivec]     ;
    h_cld_frac_i(icol_i,ipack)[ivec]      = h_cld_frac_i(icol,ipack)[ivec]     ;
    h_qv_prev(icol_i,ipack)[ivec]         = h_qv_prev(icol,ipack)[ivec]        ;
    h_T_prev(icol_i,ipack)[ivec]          = h_T_prev(icol,ipack)[ivec]         ;
    }
  }

//  Kokkos::parallel_for(
//    "iter columns",
//    ncol,
//    KOKKOS_LAMBDA(const int &i) {
//    // For column i = 0, use the ICs as originally coded in python and
//    // subsequently modified here. For columns i > 0, introduce some small
//    // variations.
//    for (int k = 0; k < nk; ++k) 
//    {
//      int ipack = k / Spack::n;
//      int ivec  = k % Spack::n;
//      int ipack_m20 = (nk-20+k) / Spack::n;
//      int ivec_m20  = (nk-20+k) % Spack::n;
//    // TODO : AaronDonahue: you were going to gather all of the k for loops together and 
//    //        use the exampler from atmosphere_microphysics.cpp to assign values using ipack and ilev, instead of i,k.
//    // max cld at ~700mb, decreasing to 0 at 900mb.
//      if (k < 15) 
//      {
//        //device_views.at("qc")(i,ipack_m20)[ivec_m20] = 1e-4*(1 - double(k)/14);
//        //device_views.at("qr")(i,ipack_m20)[ivec_m20] = 1e-5*(1 - double(k)/19);
//        d_qc(i,ipack_m20)[ivec_m20] = 1e-4*(1 - double(k)/14);
//      }
//      if (k < 20)
//      {
//        d_qr(i,ipack_m20)[ivec_m20] = 1e-5*(1 - double(k)/19);
//      }
//      d_nc(i,ipack)[ivec] = 1e6;
//      d_nr(i,ipack)[ivec] = 1e6;
//    //                                                      v (in the python)
//      if (k < 15)
//      {
//        d_qi(i,ipack_m20)[ivec_m20] = 1e-4; //*(1 - double(k)/14)
//        d_qm(i,ipack_m20)[ivec_m20] = 1e-4*(1 - double(k)/14);
//        d_bm(i,ipack_m20)[ivec_m20] = 1e-2;
//      }
//      d_ni(i,ipack)[ivec] = 1e6;
//    // guess at reasonable value based on: m3/kg is 1/density and liquid water has
//    // a density of 1000 kg/m3
//
//    // qv goes to zero halfway through profile (to avoid condensate near model
//    // top)
//      const auto tmp = -5e-4 + 1e-3/double(nk)*k;
//      d_qv(i,ipack)[ivec] = tmp > 0 ? tmp : 0;
//    // pres is actually an input variable, but needed here to compute theta.
//      d_pmid(i,ipack)[ivec] = double(100) + 1e5/double(nk)*k;
//      d_dp(i,ipack)[ivec] = 1e5/double(nk);
//      d_exner(i,ipack)[ivec] = std::pow((1e5/d_pmid(i,ipack)[ivec]), (287.15/1005.0));
//      d_cld_frac_i(i,ipack)[ivec] = 1.0;
//      d_cld_frac_l(i,ipack)[ivec] = 1.0;
//      d_cld_frac_r(i,ipack)[ivec] = 1.0;
//    // inv_qc_relvar=mean(qc)/var(qc) measures subgrid qc variability. It is computed in SHOC
//    // and used by P3. It can range between 0.1 and 10.0. Setting to a typical value of 1.0
//    // here.
//      d_inv_qc_relvar(i,ipack)[ivec] =  1.0;
//
//    // To get potential temperature, start by making absolute temperature vary
//    // between 150K at top of atmos and 300k at surface, then convert to potential
//    // temp.
//      d_T_atm(i,ipack)[ivec] = 150 + 150/double(nk)*k;
//      d_th_atm(i,ipack)[ivec] = d_T_atm(i,ipack)[ivec]*std::pow(Real(consts::P0/d_pmid(i,ipack)[ivec]), Real(consts::RD/consts::CP)); //TODO: Delete, should be handled locally
//    }
//    for (int k=0;k < 15;k++) {
//      int ipack_m20 = (nk-20+k) / Spack::n;
//      int ivec_m20  = (nk-20+k) % Spack::n;
//      d_qv(i,ipack_m20)[ivec_m20] = 5e-3;
//    }
//    int ipack_nkm1 = (nk-1) / Spack::n;
//    int ivec_nkm1  = (nk-1) % Spack::n;
//    int ipack_nkm2 = (nk-2) / Spack::n;
//    int ivec_nkm2  = (nk-2) % Spack::n;
//    int ipack, ivec;
//    // The next section modifies inout variables to satisfy weird conditions
//    // needed for code coverage.
//    d_qi(i,ipack_nkm1)[ivec_nkm1] = 1e-9;
//    d_qv(i,ipack_nkm1)[ivec_nkm1] = 5e-2; // also needs to be supersaturated to avoid getting set
//    // to 0 earlier.
//
//    // make lowest-level qc and qr>0 to trigger surface rain and drizzle
//    // calculation.
//    d_qr(i,ipack_nkm1)[ivec_nkm1] = 1e-6;
//    d_qc(i,ipack_nkm1)[ivec_nkm1] = 1e-6;
//
//    // make qi>1e-8 where qr=0 to test rain collection conditional.
//    ipack = (nk-25)/Spack::n;
//    ivec  = (nk-25)%Spack::n;
//    d_qi(i,ipack)[ivec] = 5e-8;
//
//    // make qc>0 and qr>0 where T<233.15 to test homogeneous freezing.
//    ipack = 35/Spack::n;
//    ivec  = 35%Spack::n;
//    d_qc(i,ipack)[ivec] = 1e-7;
//    d_qv(i,ipack)[ivec] = 1e-6;
//
//    // deposition/condensation-freezing needs t<258.15 and >5% supersat.
//    ipack = 33/Spack::n;
//    ivec  = 33%Spack::n;
//    d_qv(i,ipack)[ivec] = 1e-4;
//
//    // compute vertical grid spacing dz (in m) from pres and theta.
//    static constexpr Real
//      g = 9.8; // gravity, m/s^2
//    for (int k = 0; k < nk; ++k) {
//      int ipack = k / Spack::n;
//      int ivec  = k % Spack::n;
//      int ipack_m1 = (k-1) / Spack::n;
//      int ivec_m1  = (k-1) % Spack::n;
//      int ipack_p1 = (k+1) / Spack::n;
//      int ivec_p1  = (k+1) % Spack::n;
//      d_qv_prev(i,ipack)[ivec] = d_qv(i,ipack)[ivec];
//      d_T_prev(i,ipack)[ivec] =  d_T_atm(i,ipack)[ivec];
//      Real plo, phi; // pressure at cell edges, Pa
//      int ivec_p = 1 / Spack::n;
//      int ipack_p = 1 % Spack::n;
//      Real tmp_p = d_pmid(i,0)[0] - 0.5*(d_pmid(i,ivec_p)[ipack_p] - d_pmid(i,0)[0])/(1 - 0);
//      if (double(i)>tmp_p) { tmp_p = double(i); }
//      plo = (k == 0  ) ?
//        tmp_p :
//        0.5*(d_pmid(i,ipack_m1)[ivec_m1] + d_pmid(i,ipack)[ivec]);
//      phi = (k == nk-1) ?
//        d_pmid(i,ipack_nkm1)[ivec_nkm1] + 0.5*(d_pmid(i,ipack_nkm1)[ivec_nkm1] - d_pmid(i,ipack_nkm2)[ivec_nkm2])/(1 - 0) :
//        0.5*(d_pmid(i,ipack)[ivec] + d_pmid(i,ipack_p1)[ivec_p1]);
//      const auto dp = phi - plo;
//      d_zi(i,ipack)[ivec] = log(consts::P0/d_pmid(i,ipack)[ivec])*consts::RD*d_T_atm(i,ipack)[ivec]/g/0.029;
//      d_dz(i,ipack)[ivec] = consts::RD*d_T_atm(i,ipack)[ivec]/(g*d_pmid(i,ipack)[ivec])*dp; //TODO: Delete, should be handled locally
//    }
//  }); // iter column Kokkos Loop

  // Copy Host Views back to Device
  Kokkos::deep_copy(d_T_atm          , h_T_atm          ); 
  Kokkos::deep_copy(d_ast            , h_ast            ); 
  Kokkos::deep_copy(d_ni_activated   , h_ni_activated   ); 
  Kokkos::deep_copy(d_nc_nuceat_tend , h_nc_nuceat_tend ); 
  Kokkos::deep_copy(d_pmid           , h_pmid           ); 
  Kokkos::deep_copy(d_dp             , h_dp             ); 
  Kokkos::deep_copy(d_zi             , h_zi             ); 
  Kokkos::deep_copy(d_qv_prev        , h_qv_prev        ); 
  Kokkos::deep_copy(d_T_prev         , h_T_prev         ); 
  Kokkos::deep_copy(d_qv             , h_qv             ); 
  Kokkos::deep_copy(d_qc             , h_qc             ); 
  Kokkos::deep_copy(d_qr             , h_qr             ); 
  Kokkos::deep_copy(d_qi             , h_qi             ); 
  Kokkos::deep_copy(d_qm             , h_qm             ); 
  Kokkos::deep_copy(d_nc             , h_nc             ); 
  Kokkos::deep_copy(d_nr             , h_nr             ); 
  Kokkos::deep_copy(d_ni             , h_ni             ); 
  Kokkos::deep_copy(d_bm             , h_bm             ); 
  Kokkos::deep_copy(d_nccn_prescribed, h_nccn_prescribed); 
  Kokkos::deep_copy(d_inv_qc_relvar  , h_inv_qc_relvar  ); 
  // TODO: Delete eventually, should instead be set in run interface: 
  Kokkos::deep_copy(d_th_atm         , h_th_atm    ); 
  Kokkos::deep_copy(d_dz             , h_dz        ); 
  Kokkos::deep_copy(d_exner          , h_exner     ); 
  Kokkos::deep_copy(d_cld_frac_l     , h_cld_frac_l); 
  Kokkos::deep_copy(d_cld_frac_i     , h_cld_frac_i); 
  Kokkos::deep_copy(d_cld_frac_r     , h_cld_frac_r); 

  if (m_remapper) {
    m_remapper->registration_ends();

    m_remapper->remap(true);

    // Now we can destroy the remapper
    m_remapper = nullptr;
  }
}

} // namespace scream
