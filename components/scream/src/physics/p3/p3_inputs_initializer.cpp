#include "physics/p3/p3_inputs_initializer.hpp"

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
  using value_type           = Real;
  using device_type          = DefaultDevice;
  using view_type            = typename KokkosTypes<device_type>::template view<value_type*>;
  using host_view_type       = typename KokkosTypes<HostDevice>::template view<value_type*>;
  std::map<std::string,view_type> device_views;
  std::map<std::string,host_view_type> host_mirrors;
  for (auto name : fields_to_init)
  {
    // Get device views
    auto d_view =  m_fields.at(name).get_view();
    device_views.emplace(name,d_view);
    // Create host mirrors
    host_mirrors.emplace(name,Kokkos::create_mirror_view(d_view));
  }

  // Initialize all variables on using the host mirrors:
  
  for (int ii=0;ii<host_mirrors.at("T_atm").size();++ii)
  {
    host_mirrors.at("T_atm").data()[ii] = 273.0 + ii/1000.;
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
