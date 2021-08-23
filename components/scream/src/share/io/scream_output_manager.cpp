#include "scream_output_manager.hpp"
#include <cmath>
#include <memory>
#include "ekat/ekat_parameter_list.hpp"
#include "ekat/mpi/ekat_comm.hpp"

namespace scream
{

void OutputManager::
setup (const ekat::Comm& io_comm, const ekat::ParameterList& params,
       const std::shared_ptr<const fm_type>& field_mgr,
       const std::shared_ptr<const gm_type>& grids_mgr,
       const bool runtype_restart)
{
  m_io_comm   = io_comm;
  m_params    = params;
  m_field_mgr = field_mgr;
  m_grids_mgr = grids_mgr;
  m_runtype_restart = runtype_restart;

  // Check for model restart output
  // Restarts are a special kind of output that have unique characteristics.
  // The following constructs the parameter list to control the restart output stream.
  // The user can trigger restart output at runtime by adding the sublist
  // "Model Restart" to the SCREAM control YAML file.
  // A typical restart control entry may look like:
  // ----
  // Model Restart:
  //   Output:
  //     Frequency:       INT
  //     Frequency Units: STRING
  // ----
  // where Frequency Units is the units of the output frequency (ex. Steps, Months, etc.),
  // and Frequency>0 is the actual frequency. E.g., Frequency=2 meanse every other ${Frequency Units}

  // The output history should have the same restart frequency as the model restart output.
  int checkpoint_freq = 0;
  std::string checkpoint_freq_units = "None";
  if (m_params.isSublist("Model Restart")) {
    // Get restart parameters, and create a param list for the model restart output
    auto& out_params = m_params.sublist("Model Restart");
    make_restart_param_list(out_params);
    add_output_stream(out_params,true);

    // Gather restart frequency info, since the checkpointing info must match them
    checkpoint_freq       = out_params.sublist("Output").get<Int>("Frequency");
    checkpoint_freq_units = out_params.sublist("Output").get<std::string>("Frequency Units");
  }

  // Construct and store an output stream instance for each output request.
  // Typical output is controlled by parameter list which is stored in individual YAML files.
  // A list of those files is passed to the output manager via the parameter list.
  // In particular, check the 'Output YAML Files' sublist, and look for the entry
  // with the same name as the grid.
  // See srd/share/io/tests for examples.
  std::vector<std::string> empty;
  const auto& grid_name = m_field_mgr->get_grid()->name();
  auto& list_of_files = m_params.sublist("Output YAML Files").get(grid_name,empty);
  for (auto& it : list_of_files) {
    ekat::ParameterList out_params(it);
    parse_yaml_file(it,out_params);
    auto& checkpointing_params = out_params.sublist("Checkpointing");
    EKAT_REQUIRE_MSG (checkpointing_params.get("Frequency",checkpoint_freq)==checkpoint_freq,
        "Error! Output checkpointing frequency is different from the model restart output one.\n"
        "       In case model restart output is generated, checkpointing must have the same frequency.\n");
    EKAT_REQUIRE_MSG (checkpointing_params.get("Frequency Units",checkpoint_freq_units)==checkpoint_freq_units,
        "Error! Output checkpointing frequency units are different from the model restart output one.\n"
        "       In case model restart output is generated, checkpointing must have the same units.\n");
    add_output_stream(out_params,false);
  }
}
/*===============================================================================================*/
void OutputManager::run(util::TimeStamp& current_ts)
{
  for (size_t i=0; i<m_output_streams.size(); ++i) {
    // If the i-th stream is on a non-unique grid, we need to first run the remapper
    if (m_remappers[i]) {
      m_remappers[i]->remap(true);
    }
    m_output_streams[i]->run(current_ts);
  }
}
/*===============================================================================================*/
void OutputManager::finalize()
{
  for (auto& it : m_output_streams) {
    it->finalize();
  }
  m_output_streams.clear();
  m_io_comm = ekat::Comm();
  m_params = ekat::ParameterList("");
  m_field_mgr = nullptr;
  m_runtype_restart = false;
}
/*===============================================================================================*/
void OutputManager::
add_output_stream(const ekat::ParameterList& params, const bool model_restart_output)
{
  std::shared_ptr<const FieldManager<Real>> field_mgr;
  if (not m_field_mgr->get_grid()->is_unique()) {
    // We have to create a remapper from this grid to the unique one, so that we can do I/O.
    // On a uniquely distributed set of points.

    // Create a remapper from this fm's grid to its unique grid
    auto grid = m_field_mgr->get_grid();
    auto unique_grid = grid->get_unique_grid();
    auto remapper = m_grids_mgr->create_remapper(grid,unique_grid);

    // Register fields to remap in the remapper
    const auto& fnames = params.get<std::vector<std::string>>("Fields");
    remapper->registration_begins();
    for (const auto& name : fnames) {
      auto f_src = m_field_mgr->get_field(name);
      const auto& fid_src = f_src.get_header().get_identifier();
      remapper->register_field_from_src(fid_src);
    }
    remapper->registration_ends();

    // Create the remapped fields, that is, the fields on the unique grid,
    // by scanning the list of tgt fields of the remapper
    auto new_fm = std::make_shared<FieldManager<Real>>(m_field_mgr->get_grid()->get_unique_grid());
    new_fm->registration_begins();
    for (int i=0; i<remapper->get_num_registered_fields(); ++i) {
      // Create the field id in the tgt grid
      const auto& fid_src = remapper->get_src_field_id(i);
      auto fid_tgt = remapper->create_tgt_fid(fid_src);

      // Find out how large of a pack size we can use for the tgt field.
      // We can't find the pack size used for the src field, but we can
      // find the largest pack size that the src field allows, by taking
      // the log2 of the the last extent dim of the src field alloc props.
      const auto& f_src = m_field_mgr->get_field(fid_src.name());
      const int src_last_extent = f_src.get_header().get_alloc_properties().get_last_extent();
      const int pack_size = std::floor(std::log2(src_last_extent));
      FieldRequest req(fid_src,pack_size);
      new_fm->register_field(req);
    }
    new_fm->registration_ends();

    // The field mgr is now complete
    field_mgr = new_fm;

    // Now that we have both copies of each field (on grid and unique_grid),
    // we can bind fields in the remapper.
    for (const auto& name : fnames) {
      auto src = m_field_mgr->get_field(name);
      auto tgt = m_field_mgr->get_field(name);
      remapper->bind_field(src,tgt);
    }

    // The remapper is now complete
    m_remappers.push_back(remapper);
  } else {
    // No remapping needed
    field_mgr = m_field_mgr;
    m_remappers.push_back(nullptr);
  }
  auto output = std::make_shared<output_type>(m_io_comm,params,field_mgr,
                                              m_runtype_restart, model_restart_output);
  m_output_streams.push_back(output);

}
/*===============================================================================================*/
void OutputManager::make_restart_param_list(ekat::ParameterList& params)
{
  // Given the unique nature of restart files, this routine sets up the specific requirements.

  //TODO change this to some sort of case_name, TODO: control so suffix is .r.nc
  params.set<std::string>("Casename", "scorpio_restart_test");

  // Parse the parameters that controls this output instance.
  // Model restart themselves are instant snapshots (no averaging of any kind).
  params.set<std::string>("Averaging Type","Instant");

  params.set<Int>("Max Snapshots Per File",1);

  // Grab the restart fields from the field manager.
  // If a developer wants a field to be stored in the restart file than they must add the "restart" group tag
  // to that field when registering the field with the field manager.
  auto restart_fields = m_field_mgr->get_field_group("restart");

  // WARNING: the vector restart_fields.m_info->m_fields_names contains CaseInsensitiveString's.
  //          This can be a problem, since virtually all places will try to extract a
  //          std::vector<std::string> from the parameter list, causing a any_cast error.
  //          To fix this, create a std::vector<std::string> on the fly
  std::vector<std::string> fields_names;
  for (const auto& fn : restart_fields.m_info->m_fields_names) {
    fields_names.push_back(fn);
  }
  params.set("Fields",fields_names);
}
/*===============================================================================================*/
} // namespace scream
