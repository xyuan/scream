#include "physics/p3/p3_inputs_initializer.hpp"
#include "physics/share/physics_constants.hpp"
#include "physics/p3/p3_main_impl.hpp"
#include "share/io/scorpio_input.hpp"
#include "share/io/scream_scorpio_interface.hpp"

#include <array>

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

//  // Create device views
//  auto d_T_atm           = m_fields.at("T_atm").get_reshaped_view<Real**>();
//  auto d_ast             = m_fields.at("ast").get_reshaped_view<Real**>();
//  auto d_ni_activated    = m_fields.at("ni_activated").get_reshaped_view<Real**>();
//  auto d_nc_nuceat_tend  = m_fields.at("nc_nuceat_tend").get_reshaped_view<Real**>();
//  auto d_pmid            = m_fields.at("pmid").get_reshaped_view<Real**>();
//  auto d_dp              = m_fields.at("dp").get_reshaped_view<Real**>();
//  auto d_zi              = m_fields.at("zi").get_reshaped_view<Real**>();
//  auto d_qv_prev         = m_fields.at("qv_prev").get_reshaped_view<Real**>();
//  auto d_T_prev          = m_fields.at("T_prev").get_reshaped_view<Real**>();
//  auto d_qv              = m_fields.at("qv").get_reshaped_view<Real**>();
//  auto d_qc              = m_fields.at("qc").get_reshaped_view<Real**>();
//  auto d_qr              = m_fields.at("qr").get_reshaped_view<Real**>();
//  auto d_qi              = m_fields.at("qi").get_reshaped_view<Real**>();
//  auto d_qm              = m_fields.at("qm").get_reshaped_view<Real**>();
//  auto d_nc              = m_fields.at("nc").get_reshaped_view<Real**>();
//  auto d_nr              = m_fields.at("nr").get_reshaped_view<Real**>();
//  auto d_ni              = m_fields.at("ni").get_reshaped_view<Real**>();
//  auto d_bm              = m_fields.at("bm").get_reshaped_view<Real**>();
//  auto d_nccn_prescribed = m_fields.at("nccn_prescribed").get_reshaped_view<Real**>();
//  auto d_inv_qc_relvar   = m_fields.at("inv_qc_relvar").get_reshaped_view<Real**>();
//  // TODO: Delete eventually, should instead be set in run interface: 
//  auto d_th_atm          = m_fields.at("th_atm").get_reshaped_view<Real**>();  
//  auto d_dz              = m_fields.at("dz").get_reshaped_view<Real**>();  
//  auto d_exner           = m_fields.at("exner").get_reshaped_view<Real**>();  
//  auto d_cld_frac_l      = m_fields.at("cld_frac_l").get_reshaped_view<Real**>();  
//  auto d_cld_frac_i      = m_fields.at("cld_frac_i").get_reshaped_view<Real**>();  
//  auto d_cld_frac_r      = m_fields.at("cld_frac_r").get_reshaped_view<Real**>(); 
//  // Create host mirrors 
//  auto h_T_atm           = Kokkos::create_mirror_view(d_T_atm);
//  auto h_ast             = Kokkos::create_mirror_view(d_ast);
//  auto h_ni_activated    = Kokkos::create_mirror_view(d_ni_activated);
//  auto h_nc_nuceat_tend  = Kokkos::create_mirror_view(d_nc_nuceat_tend);
//  auto h_pmid            = Kokkos::create_mirror_view(d_pmid);
//  auto h_dp              = Kokkos::create_mirror_view(d_dp);
//  auto h_zi              = Kokkos::create_mirror_view(d_zi);
//  auto h_qv_prev         = Kokkos::create_mirror_view(d_qv_prev);
//  auto h_T_prev          = Kokkos::create_mirror_view(d_T_prev);
//  auto h_qv              = Kokkos::create_mirror_view(d_qv);
//  auto h_qc              = Kokkos::create_mirror_view(d_qc);
//  auto h_qr              = Kokkos::create_mirror_view(d_qr);
//  auto h_qi              = Kokkos::create_mirror_view(d_qi);
//  auto h_qm              = Kokkos::create_mirror_view(d_qm);
//  auto h_nc              = Kokkos::create_mirror_view(d_nc);
//  auto h_nr              = Kokkos::create_mirror_view(d_nr);
//  auto h_ni              = Kokkos::create_mirror_view(d_ni);
//  auto h_bm              = Kokkos::create_mirror_view(d_bm);
//  auto h_nccn_prescribed = Kokkos::create_mirror_view(d_nccn_prescribed);
//  auto h_inv_qc_relvar   = Kokkos::create_mirror_view(d_inv_qc_relvar);
//  // TODO: Delete eventually, should instead be set in run interface: 
//  auto h_th_atm          = Kokkos::create_mirror_view(d_th_atm);  
//  auto h_dz              = Kokkos::create_mirror_view(d_dz);  
//  auto h_exner           = Kokkos::create_mirror_view(d_exner);  
//  auto h_cld_frac_l      = Kokkos::create_mirror_view(d_cld_frac_l);  
//  auto h_cld_frac_i      = Kokkos::create_mirror_view(d_cld_frac_i);  
//  auto h_cld_frac_r      = Kokkos::create_mirror_view(d_cld_frac_r);  
  // Initialize all variables on using the host mirrors:
  // For now use dummy values copied from `p3_ic_cases.cpp`
  using consts          = scream::physics::Constants<Real>;
  auto temp = m_fields.at("qc");
  auto mdims = temp.get_header().get_identifier().get_layout();
  Int ncol = mdims.dim(0); 
  Int nk   = mdims.dim(1);
//  Kokkos::parallel_for(
//    "init p3 vals",
//    1,
//    KOKKOS_LAMBDA(const int& i_dum) {

//ASD  Kokkos::parallel_for(
//ASD    "iter columns",
//ASD    ncol,
//ASD    KOKKOS_LAMBDA(const int &i) {
//ASD//ASD  for (Int i = 0; i < ncol; ++i) {
//ASD    // For column i = 0, use the ICs as originally coded in python and
//ASD    // subsequently modified here. For columns i > 0, introduce some small
//ASD    // variations.
//ASD    for (int k = 0; k < nk; ++k) 
//ASD    {
//ASD      int ipack = k / Spack::n;
//ASD      int ivec  = k % Spack::n;
//ASD      int ipack_m20 = (nk-20+k) / Spack::n;
//ASD      int ivec_m20  = (nk-20+k) % Spack::n;
//ASD    // TODO : AaronDonahue: you were going to gather all of the k for loops together and 
//ASD    //        use the exampler from atmosphere_microphysics.cpp to assign values using ipack and ilev, instead of i,k.
//ASD    // max cld at ~700mb, decreasing to 0 at 900mb.
//ASD      if (k < 15) 
//ASD      {
//ASD        host_mirrors.at("qc")(i,ipack_m20)[ivec_m20] = 1e-4*(1 - double(k)/14);
//ASD        host_mirrors.at("qr")(i,ipack_m20)[ivec_m20] = 1e-5*(1 - double(k)/19);
//ASD      }
//ASD      host_mirrors.at("nc")(i,ipack)[ivec] = 1e6;
//ASD      host_mirrors.at("nr")(i,ipack)[ivec] = 1e6;
//ASD    //                                                      v (in the python)
//ASD      if (k < 15)
//ASD      {
//ASD        host_mirrors.at("qi")(i,ipack_m20)[ivec_m20] = 1e-4; //*(1 - double(k)/14)
//ASD        host_mirrors.at("qm")(i,ipack_m20)[ivec_m20] = 1e-4*(1 - double(k)/14);
//ASD        host_mirrors.at("bm")(i,ipack_m20)[ivec_m20] = 1e-2;
//ASD      }
//ASD      host_mirrors.at("ni")(i,ipack)[ivec] = 1e6;
//ASD    // guess at reasonable value based on: m3/kg is 1/density and liquid water has
//ASD    // a density of 1000 kg/m3
//ASD
//ASD    // qv goes to zero halfway through profile (to avoid condensate near model
//ASD    // top)
//ASD      const auto tmp = -5e-4 + 1e-3/double(nk)*k;
//ASD      host_mirrors.at("qv")(i,ipack)[ivec] = tmp > 0 ? tmp : 0;
//ASD      if (k < 15) host_mirrors.at("qv")(i,ipack_m20)[ivec_m20] = 5e-3;
//ASD    // pres is actually an input variable, but needed here to compute theta.
//ASD      host_mirrors.at("pmid")(i,ipack)[ivec] = double(100) + 1e5/double(nk)*double(k);
//ASD      host_mirrors.at("dp")(i,ipack)[ivec] = 1e5/double(nk);
//ASD      host_mirrors.at("exner")(i,ipack)[ivec] = std::pow((1e5/host_mirrors.at("pmid")(i,ipack)[ivec]), (287.15/1005.0));
//ASD      host_mirrors.at("cld_frac_i")(i,ipack)[ivec] = 1.0;
//ASD      host_mirrors.at("cld_frac_l")(i,ipack)[ivec] = 1.0;
//ASD      host_mirrors.at("cld_frac_r")(i,ipack)[ivec] = 1.0;
//ASD    // inv_qc_relvar=mean(qc)/var(qc) measures subgrid qc variability. It is computed in SHOC
//ASD    // and used by P3. It can range between 0.1 and 10.0. Setting to a typical value of 1.0
//ASD    // here.
//ASD      host_mirrors.at("inv_qc_relvar")(i,ipack)[ivec] = double(i)+double(k)/100.0;
//ASD
//ASD    // To get potential temperature, start by making absolute temperature vary
//ASD    // between 150K at top of atmos and 300k at surface, then convert to potential
//ASD    // temp.
//ASD      host_mirrors.at("T_atm")(i,ipack)[ivec] = 150 + 150/double(nk)*double(k);
//ASD//      if (i > 0) host_mirrors.at("T_atm")(i,ipack)[ivec] += ((i % 3) - 0.5)/double(nk)*k;
//ASD      auto pmid_tmp = host_mirrors.at("pmid")(i,ipack)[ivec];
//ASD      auto tatm_tmp = host_mirrors.at("T_atm")(i,ipack)[ivec];
//ASD      auto thatm_tmp = host_mirrors.at("th_atm")(i,ipack)[ivec];
//ASD      auto p1_tmp = Real(consts::P0/host_mirrors.at("pmid")(i,ipack)[ivec]);
//ASD      auto p2_tmp = Real(consts::RD/consts::CP);
//ASD      printf("ASD [%d] - (%2d,%2d) -> (%2d,%2d) and (%2d,%2d) :: %f, %e, %e, %e, %e, %e\n",Spack::n, 
//ASD         i,k,
//ASD         ipack,ivec,
//ASD         ipack_m20,ivec_m20,
//ASD         host_mirrors.at("inv_qc_relvar")(i,ipack)[ivec],
//ASD         pmid_tmp,
//ASD         tatm_tmp,
//ASD         thatm_tmp,
//ASD         p1_tmp,
//ASD         p2_tmp
//ASD         );
//ASD      host_mirrors.at("th_atm")(i,ipack)[ivec] = host_mirrors.at("T_atm")(i,ipack)[ivec]*std::pow(Real(consts::P0/host_mirrors.at("pmid")(i,ipack)[ivec]), Real(consts::RD/consts::CP)); //TODO: Delete, should be handled locally
//ASD    }
//ASD    int ipack_nkm1 = (nk-1) / Spack::n;
//ASD    int ivec_nkm1  = (nk-1) % Spack::n;
//ASD    int ipack_nkm2 = (nk-2) / Spack::n;
//ASD    int ivec_nkm2  = (nk-2) % Spack::n;
//ASD    // The next section modifies inout variables to satisfy weird conditions
//ASD    // needed for code coverage.
//ASD    host_mirrors.at("qi")(i,ipack_nkm1)[ivec_nkm1] = 1e-9;
//ASD    host_mirrors.at("qv")(i,ipack_nkm1)[ivec_nkm1] = 5e-2; // also needs to be supersaturated to avoid getting set
//ASD    // to 0 earlier.
//ASD
//ASD    // make lowest-level qc and qr>0 to trigger surface rain and drizzle
//ASD    // calculation.
//ASD    host_mirrors.at("qr")(i,ipack_nkm1)[ivec_nkm1] = 1e-6;
//ASD    host_mirrors.at("qc")(i,ipack_nkm1)[ivec_nkm1] = 1e-6;
//ASD
//ASD    // make qi>1e-8 where qr=0 to test rain collection conditional.
//ASD    host_mirrors.at("qi")(i,(nk-25)/Spack::n)[(nk-25)%Spack::n] = 5e-8;
//ASD
//ASD    // make qc>0 and qr>0 where T<233.15 to test homogeneous freezing.
//ASD    host_mirrors.at("qc")(i,35/Spack::n)[35%Spack::n] = 1e-7;
//ASD    host_mirrors.at("qv")(i,35/Spack::n)[35%Spack::n] = 1e-6;
//ASD
//ASD    // deposition/condensation-freezing needs t<258.15 and >5% supersat.
//ASD    host_mirrors.at("qv")(i,33/Spack::n)[33%Spack::n] = 1e-4;
//ASD
//ASD    // compute vertical grid spacing dz (in m) from pres and theta.
//ASD    static constexpr double
//ASD      g = 9.8; // gravity, m/s^2
//ASD    for (int k = 0; k < nk; ++k) {
//ASD      int ipack = k / Spack::n;
//ASD      int ivec  = k % Spack::n;
//ASD      int ipack_m1 = (k-1) / Spack::n;
//ASD      int ivec_m1  = (k-1) % Spack::n;
//ASD      int ipack_p1 = (k+1) / Spack::n;
//ASD      int ivec_p1  = (k+1) % Spack::n;
//ASD      host_mirrors.at("qv_prev")(i,ipack)[ivec] = host_mirrors.at("qv")(i,ipack)[ivec];
//ASD      host_mirrors.at("T_prev")(i,ipack)[ivec] = host_mirrors.at("T_atm")(i,ipack)[ivec];
//ASD      double plo, phi; // pressure at cell edges, Pa
//ASD      plo = (k == 0  ) ?
//ASD        std::max<double>(i, host_mirrors.at("pmid")(i,0)[0] - 0.5*(host_mirrors.at("pmid")(i,0)[1] - host_mirrors.at("pmid")(i,0)[0])/(1 - 0)) :
//ASD        0.5*(host_mirrors.at("pmid")(i,ipack_m1)[ivec_m1] + host_mirrors.at("pmid")(i,ipack)[ivec]);
//ASD      phi = (k == nk-1) ?
//ASD        host_mirrors.at("pmid")(i,ipack_nkm1)[ivec_nkm1] + 0.5*(host_mirrors.at("pmid")(i,ipack_nkm1)[ivec_nkm1] - host_mirrors.at("pmid")(i,ipack_nkm2)[ivec_nkm2])/(1 - 0) :
//ASD        0.5*(host_mirrors.at("pmid")(i,ipack)[ivec] + host_mirrors.at("pmid")(i,ipack_p1)[ivec_p1]);
//ASD      const auto dp = phi - plo;
//ASD      host_mirrors.at("zi")(i,ipack)[ivec] = log(consts::P0/host_mirrors.at("pmid")(i,ipack)[ivec])*consts::RD*host_mirrors.at("T_atm")(i,ipack)[ivec]/g/0.029;
//ASD      host_mirrors.at("dz")(i,ipack)[ivec] = consts::RD*host_mirrors.at("T_atm")(i,ipack)[ivec]/(g*host_mirrors.at("pmid")(i,ipack)[ivec])*dp; //TODO: Delete, should be handled locally
//ASD    }
//ASD//ASD  } // i for loop
//ASD  }); // iter column Kokkos Loop
  Kokkos::fence();

  // Deep copy from host view back to device view
  for (auto name : fields_to_init)
  {
    Kokkos::deep_copy(device_views.at(name),host_mirrors.at(name));
  }

  if (m_remapper) {
    m_remapper->registration_ends();

    m_remapper->remap(true);

    // Now we can destroy the remapper
    m_remapper = nullptr;
  }
}

} // namespace scream
