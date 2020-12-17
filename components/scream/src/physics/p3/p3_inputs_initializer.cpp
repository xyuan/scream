#include "physics/p3/p3_inputs_initializer.hpp"
#include "physics/share/physics_constants.hpp"

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
  using value_type      = Real;
  using device_type     = DefaultDevice;
  using DeviceView      = ekat::Unmanaged<typename KokkosTypes<DefaultDevice>::template view<value_type**> >;
  using HostView        = ekat::Unmanaged<typename KokkosTypes<HostDevice>::template view<value_type**> >;

  std::map<std::string,DeviceView> device_views;
  std::map<std::string,HostView> host_mirrors;
  for (auto name : fields_to_init)
  {
    // Get device views
    DeviceView d_view =  m_fields.at(name).get_reshaped_view<Real**>();
    device_views.emplace(name,d_view);
    // Create host mirrors
    host_mirrors.emplace(name,Kokkos::create_mirror_view(d_view));
  }

  // Initialize all variables on using the host mirrors:
  // For now use dummy values copied from `p3_ic_cases.cpp`
  using consts          = scream::physics::Constants<Real>;
  auto temp = m_fields.at("qc");
  auto mdims = temp.get_header().get_identifier().get_layout();
  Int ncol = mdims.dim(0); 
  Int nk   = mdims.dim(1); 
  for (Int i = 0; i < ncol; ++i) {
    // For column i = 0, use the ICs as originally coded in python and
    // subsequently modified here. For columns i > 0, introduce some small
    // variations.
    int k;
    // max cld at ~700mb, decreasing to 0 at 900mb.
    for (k = 0; k < 15; ++k) host_mirrors.at("qc")(i,nk-20+k) = 1e-4*(1 - double(k)/14);
    for (k = 0; k < nk; ++k) host_mirrors.at("nc")(i,k) = 1e6;
    // max rain at 700mb, decreasing to zero at surf.
    for (k = 0; k < 20; ++k) host_mirrors.at("qr")(i,nk-20+k) = 1e-5*(1 - double(k)/19);
    for (k = 0; k < nk; ++k) host_mirrors.at("nr")(i,k) = 1e6;

    //                                                      v (in the python)
    for (k = 0; k < 15; ++k) host_mirrors.at("qi")(i,nk-20+k) = 1e-4; //*(1 - double(k)/14)
    for (k = 0; k < nk; ++k) host_mirrors.at("ni")(i,k) = 1e6;
    for (k = 0; k < 15; ++k) host_mirrors.at("qm")(i,nk-20+k) = 1e-4*(1 - double(k)/14);
    // guess at reasonable value based on: m3/kg is 1/density and liquid water has
    // a density of 1000 kg/m3
    for (k = 0; k < 15; ++k) host_mirrors.at("bm")(i,nk-20+k) = 1e-2;

    // qv goes to zero halfway through profile (to avoid condensate near model
    // top)
    for (k = 0; k < nk; ++k) {
      const auto tmp = -5e-4 + 1e-3/double(nk)*k;
      host_mirrors.at("qv")(i,k) = tmp > 0 ? tmp : 0;
    }
    // make layer with qc saturatehost_mirrors.at("")
    for (k = 0; k < 15; ++k) host_mirrors.at("qv")(i,nk-20+k) = 5e-3;

    // pres is actually an input variable, but needed here to compute theta.
    for (k = 0; k < nk; ++k) host_mirrors.at("pmid")(i,k) = 100 + 1e5/double(nk)*k;
    // dp is actually an input variable, but needed here to compute theta.
    for (k = 0; k < nk; ++k) host_mirrors.at("dp")(i,k) = 1e5/double(nk);
    // inv_qc_relvar=mean(qc)/var(qc) measures subgrid qc variability. It is computed in SHOC
    // and used by P3. It can range between 0.1 and 10.0. Setting to a typical value of 1.0
    // here.
    for (k = 0; k < nk; ++k) host_mirrors.at("inv_qc_relvar")(i,k) = 1.;

    // To get potential temperature, start by making absolute temperature vary
    // between 150K at top of atmos and 300k at surface, then convert to potential
    // temp.
    for (k = 0; k < nk; ++k) {
      host_mirrors.at("T_atm")(i,k) = 150 + 150/double(nk)*k;
      if (i > 0) host_mirrors.at("T_atm")(i,k) += ((i % 3) - 0.5)/double(nk)*k;
    }

    // The next section modifies inout variables to satisfy weird conditions
    // needed for code coverage.
    host_mirrors.at("qi")(i,nk-1) = 1e-9;
    host_mirrors.at("qv")(i,nk-1) = 5e-2; // also needs to be supersaturated to avoid getting set
    // to 0 earlier.

    // make lowest-level qc and qr>0 to trigger surface rain and drizzle
    // calculation.
    host_mirrors.at("qr")(i,nk-1) = 1e-6;
    host_mirrors.at("qc")(i,nk-1) = 1e-6;

    // make qi>1e-8 where qr=0 to test rain collection conditional.
    host_mirrors.at("qi")(i,nk-25) = 5e-8;

    // make qc>0 and qr>0 where T<233.15 to test homogeneous freezing.
    host_mirrors.at("qc")(i,35) = 1e-7;
    host_mirrors.at("qv")(i,35) = 1e-6;

    // deposition/condensation-freezing needs t<258.15 and >5% supersat.
    host_mirrors.at("qv")(i,33) = 1e-4;

    // set qv_prev and t_prev to qv and T vals
    for (k = 0; k < nk; ++k){
      host_mirrors.at("qv_prev")(i,k) = host_mirrors.at("qv")(i,k);
      host_mirrors.at("T_prev")(i,k) = host_mirrors.at("T_atm")(i,k);
    }

    // compute vertical grid spacing dz (in m) from pres and theta.
    static constexpr double
      g = 9.8; // gravity, m/s^2
    for (k = 0; k < nk; ++k) {
      double plo, phi; // pressure at cell edges, Pa
      plo = (k == 0   ) ?
        std::max<double>(i, host_mirrors.at("pmid")(i,0) - 0.5*(host_mirrors.at("pmid")(i,1) - host_mirrors.at("pmid")(i,0))/(1 - 0)) :
        0.5*(host_mirrors.at("pmid")(i,k-1) + host_mirrors.at("pmid")(i,k));
      phi = (k == nk-1) ?
        host_mirrors.at("pmid")(i,nk-1) + 0.5*(host_mirrors.at("pmid")(i,nk-1) - host_mirrors.at("pmid")(i,nk-2))/(1 - 0) :
        0.5*(host_mirrors.at("pmid")(i,k) + host_mirrors.at("pmid")(i,k+1));
      const auto dp = phi - plo;
      host_mirrors.at("zi")(i,k)   = log(consts::P0/host_mirrors.at("pmid")(i,k))*consts::RD*host_mirrors.at("T_atm")(i,k)/g/0.029;
    }
  }

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
