#include "control/atmosphere_driver.hpp"

#include "share/grid/remap//abstract_remapper.hpp"
#include "share/atm_process/atmosphere_process.hpp"

#include "ekat/ekat_pack.hpp"

namespace scream {

// === A dummy atm process, on Point grid === //

class DummyProcess : public scream::AtmosphereProcess {
public:

  DummyProcess (const ekat::Comm& comm, const ekat::ParameterList& params)
   : m_comm(comm)
  {
    m_params = params;
  }

  // The type of the block: Physics
  AtmosphereProcessType type () const { return AtmosphereProcessType::Physics; }

  std::set<std::string> get_required_grids () const {
    return std::set<std::string> {m_params.get<std::string>("Grid Name")};
  }

  // Return some sort of name, linked to PType
  std::string name () const { return m_name; }

  // The communicator associated with this atm process
  const ekat::Comm& get_comm () const { return m_comm; }

  void set_grids (const std::shared_ptr<const GridsManager> grids_manager) {
    using namespace ShortFieldTagsNames;

    // Grid and ref grid
    m_grid = grids_manager->get_grid(m_params.get<std::string>("Grid Name"));
    auto ref_grid = grids_manager->get_reference_grid();

    // If grid==ref_grid, these are IdentityRemapper instances
    m_remap_in  = grids_manager->create_remapper_from_ref_grid(m_grid);
    m_remap_out = grids_manager->create_remapper_to_ref_grid(m_grid);
    m_remap_in->registration_begins();
    m_remap_out->registration_begins();

    // Grid specs
    m_ncols = m_grid->get_num_local_dofs();
    m_nlevs = m_grid->get_num_vertical_levels();

    // Create fids and group requests
    FieldLayout layout = ref_grid->get_3d_scalar_layout(true);
    auto units = ekat::units::m;

    const auto& fields = m_params.get<std::vector<std::string>>("Fields");
    const auto& groups = m_params.get<std::vector<std::string>>("Groups");

    for (const auto& name : groups) {
      const auto& pl = m_params.sublist(name);

      GroupRequest req(name,ref_grid->name(),SCREAM_PACK_SIZE,false);
      if (pl.isParameter("Super Group")) {
        req.superset_group = pl.get<std::string>("Super Group");
        for (const auto& n : pl.get<std::vector<std::string>>("Exclude")) {
          req.exclude_superset_fields.push_back(n);
        }
      } else {
        auto members = pl.get<std::vector<std::string>>("Members");
        for (const auto& n : members) {
          m_f2g[n].insert(name);
        }
      }
      m_groups_req.push_back(req);
    }

    for (const auto& name : fields) {
      m_fids.emplace(name,layout,units,ref_grid->name());
    }
  }

  // Register all fields in the given repo
  void register_fields (FieldRepository<Real>& field_repo) const {
    for (const auto& fid : m_fids) {
      if (m_f2g.find(fid.name())!=m_f2g.end()) {
        field_repo.register_field(fid,m_f2g.at(fid.name()));
      } else {
        field_repo.register_field(fid);
      }
    }
  }

  void set_updated_group (const FieldGroup<Real>& field_group) {
    for (const auto& it : field_group.m_fields) {
      const auto& f = it.second;
      const auto& fid = f->get_header().get_identifier();
      m_ref_fields[fid.name()] = *f;
    }
    m_groups.emplace(field_group.m_info->m_group_name,field_group);
  }

  // Providing a list of required and computed fields
  const std::set<FieldIdentifier>&  get_required_fields () const { return m_fids; }
  const std::set<FieldIdentifier>&  get_computed_fields () const { return m_fids; }

  std::list<GroupRequest> get_updated_groups () const {
    return m_groups_req;
  }

protected:

  void initialize_impl (const util::TimeStamp&) {
    // Set fields in the remappers
    for (auto& it : m_ref_fields) {
      auto ref_fid = it.second.get_header().get_identifier();
      auto fid = m_remap_in->create_tgt_fid(ref_fid);

      Field<Real> f(fid);
      f.allocate_view();
      m_fields[it.first] = f;

      m_remap_in->register_field(it.second,f);
      m_remap_out->register_field(f,it.second);
    }

    m_remap_in->registration_ends();
    m_remap_out->registration_ends();
  }

// CUDA needs top level lambdas to be enclosed by a method that is public.
#ifdef KOKKOS_ENABLE_CUDA
public:
#endif
  void run_impl (const Real /* dt */) {
    m_remap_in->remap(true);
    for (auto it : m_fields) {
      it.second.sync_to_host();
      auto vh = it.second.get_reshaped_view<Real**,Host>();
      if (m_grid_name=="Col Lev") {
        for (int i=0; i<m_ncols; ++i) {
          for (int k=0; k<m_nlevs; ++k) {
            vh(i,k) *= 2.0;
          }
        }
      }
      it.second.sync_to_dev();
    }
    m_remap_out->remap(true);
  }

protected:

  void finalize_impl () {
    // Do nothing
  }

  // Setting the field in the atmosphere process
  void set_required_field_impl (const Field<const Real>& /* f */) {
    // All fields are in/out, so wait till 'set_computed_field_impl',
    // so we can store as non const field
  }
  void set_computed_field_impl (const Field<Real>& f) {
    m_ref_fields[f.get_header().get_identifier().name()] = f;
  }

  std::set<FieldIdentifier> m_fids;
  std::list<GroupRequest>   m_groups_req;

  std::map<std::string,Field<Real>>       m_ref_fields;
  std::map<std::string,Field<Real>>       m_fields;
  std::map<std::string,FieldGroup<Real>>  m_groups;

  std::map<std::string,std::set<std::string>> m_f2g;

  std::shared_ptr<AbstractRemapper<Real>> m_remap_in;
  std::shared_ptr<AbstractRemapper<Real>> m_remap_out;

  std::shared_ptr<const AbstractGrid>   m_grid;

  std::string m_name;
  std::string m_grid_name;

  int m_ncols;
  int m_nlevs;

  ekat::ParameterList m_params;

  ekat::Comm    m_comm;
};

} // namespace scream
