#include <catch2/catch.hpp>

#include "control/atmosphere_driver.hpp"

#include "physics/p3/atmosphere_microphysics.hpp"

#include "physics/share/physics_only_grids_manager.hpp"

#include "share/atm_process/atmosphere_process.hpp"

#include "ekat/ekat_pack.hpp"
#include "ekat/ekat_parse_yaml_file.hpp"

#include "physics/p3/p3_ic_cases.hpp"
#include "physics/p3/p3_f90.hpp"
#include "physics/p3/p3_functions_f90.hpp"
namespace scream {

// === A dummy physics grids for this test === //

TEST_CASE("p3-stand-alone", "") {
  using namespace scream;
  using namespace scream::control;
  using namespace scream::p3;

  constexpr int num_iters = 1;
  constexpr Real dt=300.0;

  // Load ad parameter list
  std::string fname = "input.yaml";
  ekat::ParameterList ad_params("Atmosphere Driver");
  REQUIRE_NOTHROW ( parse_yaml_file(fname,ad_params) );

  // Create a comm
  ekat::Comm atm_comm (MPI_COMM_WORLD);

  // Need to register products in the factory *before* we create any atm process or grids manager.,
  auto& proc_factory = AtmosphereProcessFactory::instance();
  auto& gm_factory = GridsManagerFactory::instance();
  proc_factory.register_product("p3",&create_atmosphere_process<P3Microphysics>);
  gm_factory.register_product("Physics Only",&physics::create_physics_only_grids_manager);

  // Create the grids manager
  auto& gm_params = ad_params.sublist("Grids Manager");
  const std::string& gm_type = gm_params.get<std::string>("Type");
  auto gm = GridsManagerFactory::instance().create(gm_type,atm_comm,gm_params);

  // Create the driver
  AtmosphereDriver ad;

  // Init and run (do not finalize, or you'll clear the field repo!)
  util::TimeStamp time (0,0,0,0);
  ad.initialize(atm_comm,ad_params,time);

  /* ---------------------------------------------------------------
   * setup for stand-alone test against p3_run_and_cmp baseline.
   * Check for BFB
   * --------------------------------------------------------------
   */
  // Get dimension sizes from the field manager
  const auto& grid = ad.get_grids_manager()->get_grid("Physics");
  const auto& field_mgr = *ad.get_field_mgr(grid->name());
  int ncol = grid->get_num_local_dofs();
  int nlev = grid->get_num_vertical_levels();

  // Setup test case following p3_run_and_cmp
  const auto F90_data = ic::Factory::create(ic::Factory::mixed, ncol);
  F90_data->dt = dt;
  F90_data->it = num_iters;
  F90_data->do_predict_nc = false;
  F90_data->do_prescribed_CCN = false;

  p3_init();
  for (int i=0; i<num_iters; ++i) {
     p3_main(*F90_data, false);
  }

  // Resume running AD
  for (int i=0; i<num_iters; ++i) {
    ad.run(dt);
  }

  // TODO: get the field repo from the driver, and go get (one of)
  //       the output(s) of P3, to check its numerical value (if possible)

  // Finalize 
  ad.finalize();

  // If we got here, we were able to run p3
  REQUIRE(true);
}

} // empty namespace
