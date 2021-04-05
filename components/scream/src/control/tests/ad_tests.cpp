#include "dummy_atm_setup.hpp"

#include "control/atmosphere_driver.hpp"

#include "ekat/ekat_parameter_list.hpp"
#include "ekat/ekat_parse_yaml_file.hpp"
#include "ekat/std_meta/ekat_std_utils.hpp"

#include <catch2/catch.hpp>

namespace scream {

TEST_CASE ("group_requirements","[!throws]")
{
  constexpr int num_cols = 4;
  constexpr int num_vl   = 2;

  // Load ad parameter list
  std::string fname = "ad_tests.yaml";
  ekat::ParameterList ad_params("Atmosphere Driver");
  REQUIRE_NOTHROW ( parse_yaml_file(fname,ad_params) );

  // Make field groups info available to atm procs
  const auto& fg = ad_params.sublist("Field Groups");
  auto& ap = ad_params.sublist("Atmosphere Processes");
  for (int i=0; i<ap.get<int>("Number of Entries"); ++i) {
    auto& ap_i = ap.sublist(ekat::strint("Process",i));
    const auto& gnames = ap_i.get<std::vector<std::string>>("Groups");
    for (auto n : gnames) {
      ap_i.sublist(n) = fg.sublist(n);;
    }
  }

  // Create a comm
  ekat::Comm atm_comm (MPI_COMM_WORLD);

  // Setup the atm factories and grid manager
  dummy_atm_init(num_cols, num_vl, atm_comm);

  // Create the driver
  control::AtmosphereDriver ad;

  // Init and run a single time step
  util::TimeStamp init_time(0,0,0,0.0);
  ad.initialize(atm_comm,ad_params,init_time);
  ad.run(1.0);

  // Check that field group G2 is a subset of G1
  const auto& repo = ad.get_field_repo();
  const auto& gm   = ad.get_grids_manager();

  auto G1 = repo.get_field_group("G1",gm->get_reference_grid()->name());
  auto G2 = repo.get_field_group("G2",gm->get_reference_grid()->name());
  auto G1_pl = fg.sublist("G1");
  auto G2_pl = fg.sublist("G2");
  auto G1_members= G1_pl.get<std::vector<std::string>>("Members");
  auto G2_exclude= G2_pl.get<std::vector<std::string>>("Exclude");

  // Make sure we have something in the group (avoid bugs where nothing is read correctly)
  REQUIRE (G1_members.size()>0);
  REQUIRE (G1_members.size()==G1.m_info->size());
  REQUIRE (G2.m_info->size()==(G1.m_info->size()-G2_exclude.size()));

  // Check G1 contains what we expected (and nothing more)
  for (const auto& f : G1.m_info->m_fields_names) {
    REQUIRE(ekat::contains(G1_members,f));
  }

  // Check G2 is indeed a subgroup of G1
  for (const auto& f : G2.m_info->m_fields_names) {
    REQUIRE( not ekat::contains(G2_exclude,f));
  }

  // Cleanup the AD
  ad.finalize ();

  // Cleanup atm factories and grids manager
  dummy_atm_cleanup();
}

} // namespace scream
