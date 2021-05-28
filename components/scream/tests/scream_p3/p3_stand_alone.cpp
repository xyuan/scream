#include <catch2/catch.hpp>

#include "control/atmosphere_driver.hpp"

#include "physics/p3/atmosphere_microphysics.hpp"

#include "physics/share/physics_only_grids_manager.hpp"

#include "share/atm_process/atmosphere_process.hpp"

#include "ekat/ekat_pack.hpp"
#include "ekat/ekat_parse_yaml_file.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"

#include "physics/p3/p3_ic_cases.hpp"
#include "physics/p3/p3_f90.hpp"
#include "physics/p3/p3_functions_f90.hpp"
namespace scream {

// === A dummy physics grids for this test === //

using KT = KokkosTypes<DefaultDevice>;
using ExeSpace = KT::ExeSpace;
using MemberType = KT::MemberType;
using real_view_2d = typename KT::template view_2d<Real>;

void print_sum(const int ncol, const int nlev, const std::string name, const Real* x);

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

  // Set the FM variables passed to P3 through the AD to the values set by the
  // Fortran Data Iterator
  auto pmid               = field_mgr.get_field("p_mid").get_reshaped_view<Real**>();
  auto pseudo_density     = field_mgr.get_field("pseudo_density").get_reshaped_view<Real**>();
  auto T_atm              = field_mgr.get_field("T_mid").get_reshaped_view<Real**>();
  auto cld_frac_t         = field_mgr.get_field("cldfrac_tot").get_reshaped_view<Real**>();
  auto zi                 = field_mgr.get_field("z_int").get_reshaped_view<Real**>();
  auto qv                 = field_mgr.get_field("qv").get_reshaped_view<Real**>();
  auto qv_prev            = field_mgr.get_field("qv_prev_micro_step").get_reshaped_view<Real**>();
  auto t_prev             = field_mgr.get_field("T_prev_micro_step").get_reshaped_view<Real**>();
  auto qc                 = field_mgr.get_field("qc").get_reshaped_view<Real**>();
  auto nc                 = field_mgr.get_field("nc").get_reshaped_view<Real**>();
  auto qr                 = field_mgr.get_field("qr").get_reshaped_view<Real**>();
  auto nr                 = field_mgr.get_field("nr").get_reshaped_view<Real**>();
  auto qi                 = field_mgr.get_field("qi").get_reshaped_view<Real**>();
  auto qm                 = field_mgr.get_field("qm").get_reshaped_view<Real**>();
  auto ni                 = field_mgr.get_field("ni").get_reshaped_view<Real**>();
  auto bm                 = field_mgr.get_field("bm").get_reshaped_view<Real**>();
  auto nc_nuceat_tend     = field_mgr.get_field("nc_nuceat_tend").get_reshaped_view<Real**>();
  auto nccn               = field_mgr.get_field("nc_activated").get_reshaped_view<Real**>();
  auto ni_activated       = field_mgr.get_field("ni_activated").get_reshaped_view<Real**>();
  auto inv_qc_relvar      = field_mgr.get_field("inv_qc_relvar").get_reshaped_view<Real**>();
  auto inv_exner          = field_mgr.get_field("inv_exner").get_reshaped_view<Real**>(); 
  auto th_atm             = field_mgr.get_field("th_atm").get_reshaped_view<Real**>(); 
  auto cld_frac_l         = field_mgr.get_field("cld_frac_l").get_reshaped_view<Real**>(); 
  auto cld_frac_i         = field_mgr.get_field("cld_frac_i").get_reshaped_view<Real**>(); 
  auto cld_frac_r         = field_mgr.get_field("cld_frac_r").get_reshaped_view<Real**>(); 
  auto dz                 = field_mgr.get_field("dz").get_reshaped_view<Real**>(); 
  {
  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(ncol, nlev);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const int i = team.league_rank();
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nlev), [&] (const int& k) {
      pmid              (i,k)  = F90_data->pres(i,k) ;
      pseudo_density    (i,k)  = F90_data->dpres(i,k) ;
      cld_frac_t        (i,k)  = F90_data->cld_frac_l(i,k) ;
      qv                (i,k)  = F90_data->qv(i,k) ;
      qv_prev           (i,k)  = F90_data->qv_prev(i,k) ;
      t_prev            (i,k)  = F90_data->t_prev(i,k) ;
      qc                (i,k)  = F90_data->qc(i,k) ;
      nc                (i,k)  = F90_data->nc(i,k) ;
      qr                (i,k)  = F90_data->qr(i,k) ;
      nr                (i,k)  = F90_data->nr(i,k) ;
      qi                (i,k)  = F90_data->qi(i,k) ;
      qm                (i,k)  = F90_data->qm(i,k) ;
      ni                (i,k)  = F90_data->ni(i,k) ;
      bm                (i,k)  = F90_data->bm(i,k) ;
      nc_nuceat_tend    (i,k)  = F90_data->nc_nuceat_tend(i,k) ;
      nccn              (i,k)  = F90_data->nccn_prescribed(i,k) ;
      ni_activated      (i,k)  = F90_data->ni_activated(i,k) ;
      inv_qc_relvar     (i,k)  = F90_data->inv_qc_relvar(i,k) ;
      inv_exner         (i,k)  = F90_data->inv_exner(i,k) ;
      th_atm            (i,k)  = F90_data->th_atm(i,k) ;
      cld_frac_l        (i,k)  = F90_data->cld_frac_l(i,k) ;
      cld_frac_i        (i,k)  = F90_data->cld_frac_i(i,k) ;
      cld_frac_r        (i,k)  = F90_data->cld_frac_r(i,k) ;
      dz                (i,k)  = F90_data->dz(i,k) ;
      // Missing: z_int, T_mid
    });
  });

  }
  Kokkos::fence();
  print_sum(ncol,nlev,"f90 qv",F90_data->qv.data());
  print_sum(ncol,nlev,"scream qv",qv.data());
  

  // Run p3_run_and_cmp type simulation
  p3_init();
  for (int i=0; i<num_iters; ++i) {
     p3_main(*F90_data, true);
     print_sum(ncol,nlev,"f90 ",F90_data->th_atm.data());
  }

  // Resume running AD
  for (int i=0; i<num_iters; ++i) {
    ad.run(dt);
    print_sum(ncol,nlev,"scream ",th_atm.data());
  }

  // TODO: get the field repo from the driver, and go get (one of)
  //       the output(s) of P3, to check its numerical value (if possible)

  // Finalize 
  ad.finalize();

  // If we got here, we were able to run p3
  REQUIRE(true);
}

void print_sum(const int ncol, const int nlev, const std::string name, const Real* x)
{
  Real result;
  Kokkos::parallel_reduce(ncol*nlev, [&] (const int& i, Real& m_sum) {
    m_sum += x[i];
  },result);
  Kokkos::fence();

  printf("Total sum for %s: %e\n",name.c_str(),result);
}

} // empty namespace
