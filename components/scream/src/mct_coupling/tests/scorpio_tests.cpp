#include <catch2/catch.hpp>

#include "scream_config.h"
#include "ekat/scream_pack.hpp"
#include "mct_coupling/scream_scorpio_interface.hpp"
#include "ekat/mpi/scream_comm.hpp"
#include "ekat/scream_types.hpp"
#include "ekat/util/ekat_md_array.hpp"



namespace {

Real f_x(const Real x, const Real t); 
Real f_y(const Real y, const Real t); 
Real f_z(const Real z, const Real t); 
Int ind_x(const Int ii);
Int ind_y(const Int jj);
Int ind_z(const Int kk);
Int ind_t(const Int tt);
void get_dof(const std::size_t dimsize, const Int myrank, const Int numranks, Int dof_len, Int &istart, Int &istop);

TEST_CASE("scorpio_interface_output", "") {

  using namespace scream;
  using namespace scream::scorpio;
  using ekat::util::data;

  // Create the set of SCORPIO output files and their respective
  // dimensions and variables.
  int compid=0;
  MPI_Fint fcomm = MPI_Comm_c2f(MPI_COMM_WORLD);
  eam_init_pio_subsystem(fcomm,compid,true);   // Gather the initial PIO subsystem data creater by component coupler
  // Register the set of output files:
  std::string outfilename = "scorpio_output_test.nc";
  register_outfile(outfilename);
  // Register the set of dimensions per output file
  int xlen=10, ylen=5, zlen=2;
  register_dimension(outfilename,"x","horizontal distance",xlen);
  register_dimension(outfilename,"y","vertical distance",ylen);
  register_dimension(outfilename,"z","height",zlen);
  register_dimension(outfilename,"time","time",0);
  // Register the set of variables per output file
  const char* vec_time[] = {"time"};
  const char* vec_x[]    = {"x"};
  const char* vec_y[]    = {"y"};
  const char* vec_z[]    = {"z"};
  const char* vec_xt[]   = {"x","time"};
  const char* vec_xyt[]  = {"x","y","time"};
  const char* vec_xyzt[] = {"x","y","z","time"};
 
  register_variable(outfilename,"time","time",1,vec_time, PIO_REAL,"t");
  register_variable(outfilename,"x","x-direction",1,vec_x, PIO_REAL,"x-real");
  register_variable(outfilename,"y","y-direction",1,vec_y, PIO_REAL,"y-real");
  register_variable(outfilename,"z","z-direction",1,vec_z, PIO_REAL,"z-real");
  register_variable(outfilename,"data_1d","test value for 1d field",2,vec_xt, PIO_REAL,"xt-real");
  register_variable(outfilename,"data_2d","test value for 2d field",3,vec_xyt, PIO_REAL,"xyt-real");
  register_variable(outfilename,"data_3d","test value for 3d field",4,vec_xyzt, PIO_REAL,"xyzt-real");
  register_variable(outfilename,"index_1d","test value for 1d field",2,vec_xt, PIO_INT,"xt-int");
  register_variable(outfilename,"index_2d","test value for 2d field",3,vec_xyt, PIO_INT,"xyt-int");
  register_variable(outfilename,"index_3d","test value for 3d field",4,vec_xyzt, PIO_INT,"xyzt-int");
  // Finished with the initialization of variables in output file
  eam_pio_enddef(outfilename);
  // Create data to be written
  std::array<Real,10> x_data;
  std::array<Real, 5> y_data;
  std::array<Real, 2> z_data;
  std::array<Int,1> xdim = {10};
  std::array<Int,1> ydim = {5};
  std::array<Int,1> zdim = {2};
  std::array<Int,1> dimlen_1d = {10};
  std::array<Int,2> dimlen_2d = {5,10};
  std::array<Int,3> dimlen_3d = {2,5,10};
  ekat::util::md_array<Real,10>       test_data_1d;
  ekat::util::md_array<Real, 5,10>    test_data_2d;
  ekat::util::md_array<Real, 2, 5,10> test_data_3d;
  ekat::util::md_array<Int,10>        test_index_1d;
  ekat::util::md_array<Int, 5,10>     test_index_2d;
  ekat::util::md_array<Int, 2, 5,10>  test_index_3d;
  Real pi = 2*acos(0.0);

  for (decltype(x_data)::size_type ii=0;ii<x_data.size();ii++) {
    x_data[ii] = 2.0*pi/x_data.size()*(ii+1);
  }
  for (int jj=0;jj<5;jj++) {
    y_data[jj] = 4.0*pi/y_data.size()*(jj+1);
  }
  for (int kk=0;kk<2;kk++) {
    z_data[kk] = 100*(kk+1);
  }
  // Write dimension data and initial fields
  grid_write_data_array(outfilename,"x",xdim,ekat::util::data(x_data));
  grid_write_data_array(outfilename,"y",ydim,ekat::util::data(y_data));
  grid_write_data_array(outfilename,"z",zdim,ekat::util::data(z_data));
  sync_outfile(outfilename); 
  // write multiple timesteps of data for comparison:
  Real dt = 1.0;
  for (int tt=0;tt<3;tt++) {
    for (decltype(x_data)::size_type ii=0;ii<x_data.size();ii++) {
      test_data_1d[ii]  = f_x(x_data[ii],tt*dt);
      test_index_1d[ii] = ind_x(ii) + ind_t(tt);
      for (int jj=0;jj<5;jj++) {
        test_data_2d[jj][ii]  = f_x(x_data[ii],tt*dt)*f_y(y_data[jj],tt*dt);
        test_index_2d[jj][ii] = ind_y(jj) + ind_x(ii) + ind_t(tt);
        for (int kk=0;kk<2;kk++) {
          test_data_3d[kk][jj][ii]  = f_x(x_data[ii],tt*dt)*f_y(y_data[jj],tt*dt) + f_z(z_data[kk],tt*dt);
          test_index_3d[kk][jj][ii] = ind_z(kk) + ind_y(jj) + ind_x(ii) + ind_t(tt);
        } //kk
      } //jj
    } //ii
    pio_update_time(outfilename,tt*dt);
    grid_write_data_array(outfilename,"index_1d",dimlen_1d,ekat::util::data(test_index_1d));
    grid_write_data_array(outfilename,"index_2d",dimlen_2d,ekat::util::data(test_index_2d));
    grid_write_data_array(outfilename,"index_3d",dimlen_3d,ekat::util::data(test_index_3d));
    grid_write_data_array(outfilename,"data_1d",dimlen_1d,ekat::util::data(test_data_1d));
    grid_write_data_array(outfilename,"data_2d",dimlen_2d,ekat::util::data(test_data_2d));
    grid_write_data_array(outfilename,"data_3d",dimlen_3d,ekat::util::data(test_data_3d));
    sync_outfile(outfilename); 
  } //tt
  ///* Now close the output file and reopen it to check that the output is correct. */
  //eam_pio_closefile(outfilename);
  //register_infile(outfilename);
  //
  //register_variable(outfilename,"time","time",1,vec_time, PIO_REAL,"t");
  //register_variable(outfilename,"x","x-direction",1,vec_x, PIO_REAL,"x-real");
  //register_variable(outfilename,"y","y-direction",1,vec_y, PIO_REAL,"y-real");
  //register_variable(outfilename,"z","z-direction",1,vec_z, PIO_REAL,"z-real");
  //register_variable(outfilename,"data_1d","test value for 1d field",1,vec_x, PIO_REAL,"x-real");
  //register_variable(outfilename,"data_2d","test value for 2d field",2,vec_xy, PIO_REAL,"xy-real");
  //register_variable(outfilename,"data_3d","test value for 3d field",3,vec_xyz, PIO_REAL,"xyz-real");
  //register_variable(outfilename,"index_1d","test value for 1d field",1,vec_x, PIO_INT,"x-int");
  //register_variable(outfilename,"index_2d","test value for 2d field",2,vec_xy, PIO_INT,"xy-int");
  //register_variable(outfilename,"index_3d","test value for 3d field",3,vec_xyz, PIO_INT,"xyz-int");


  eam_pio_finalize();
} // TEST scorpio_interface_output
/* ================================================================================================================ */
TEST_CASE("scorpio_interface_input", "") {

  using namespace scream;
  using namespace scream::scorpio;
  using ekat::util::data;

  Real pi = 2*acos(0.0);
  Real dt = 1.0;
  int compid=0;
  MPI_Fint fcomm = MPI_Comm_c2f(MPI_COMM_WORLD);
  eam_init_pio_subsystem(fcomm,compid,true);   // Gather the initial PIO subsystem data creater by component coupler
  // Register the set of output files:
  std::string infilename = "scorpio_output_baseline.nc";
  register_infile(infilename);

  const char* vec_time[] = {"time"};
  const char* vec_x[]    = {"x"};
  const char* vec_y[]    = {"y"};
  const char* vec_z[]    = {"z"};
  const char* vec_xy[]   = {"x","y"}; 
  const char* vec_xyz[]  = {"x","y","z"};
 
  register_variable(infilename,"time","time",1,vec_time, PIO_REAL,"t");
  register_variable(infilename,"x","x-direction",1,vec_x, PIO_REAL,"x-real");
  register_variable(infilename,"y","y-direction",1,vec_y, PIO_REAL,"y-real");
  register_variable(infilename,"z","z-direction",1,vec_z, PIO_REAL,"z-real");
  register_variable(infilename,"data_1d","test value for 1d field",1,vec_x, PIO_REAL,"x-real");
  register_variable(infilename,"data_2d","test value for 2d field",2,vec_xy, PIO_REAL,"xy-real");
  register_variable(infilename,"data_3d","test value for 3d field",3,vec_xyz, PIO_REAL,"xyz-real");
  register_variable(infilename,"index_1d","test value for 1d field",1,vec_x, PIO_INT,"x-int");
  register_variable(infilename,"index_2d","test value for 2d field",2,vec_xy, PIO_INT,"xy-int");
  register_variable(infilename,"index_3d","test value for 3d field",3,vec_xyz, PIO_INT,"xyz-int");

  // Create data to be written
  std::array<Int,1> xdim = {10};
  std::array<Int,1> ydim = {5};
  std::array<Int,1> zdim = {2};
  std::array<Int,1> dimlen_1d = {10};
  std::array<Int,2> dimlen_2d = {5,10};
  std::array<Int,3> dimlen_3d = {2,5,10};
  // Data to be loaded from input file
  std::array<Real,10> x_data;
  std::array<Real, 5> y_data;
  std::array<Real, 2> z_data;
  ekat::util::md_array<Real,10>       test_data_1d;
  ekat::util::md_array<Real, 5,10>    test_data_2d;
  ekat::util::md_array<Real, 2, 5,10> test_data_3d;
  ekat::util::md_array<Int,10>        test_index_1d;
  ekat::util::md_array<Int, 5,10>     test_index_2d;
  ekat::util::md_array<Int, 2, 5,10>  test_index_3d;
  // Local arrays for BFB comparison
  std::array<Real,10> comp_x;
  std::array<Real, 5> comp_y;
  std::array<Real, 2> comp_z;
  ekat::util::md_array<Real,3 ,10>       comp_data_1d;
  ekat::util::md_array<Real,3 , 5,10>    comp_data_2d;
  ekat::util::md_array<Real,3 , 2, 5,10> comp_data_3d;
  ekat::util::md_array<Int, 3,10>        comp_index_1d;
  ekat::util::md_array<Int, 3, 5,10>     comp_index_2d;
  ekat::util::md_array<Int, 3, 2, 5,10>  comp_index_3d;
  for (int tt=0;tt<3;tt++) {
    for (decltype(comp_x)::size_type ii=0;ii<x_data.size();ii++) {
      comp_x[ii] = 2.0*pi/comp_x.size()*(ii+1);
      comp_data_1d[tt][ii] = 0.1 * cos(comp_x[ii]+tt*dt);
      comp_index_1d[tt][ii] = ii + 10000*tt;
      for (decltype(comp_y)::size_type jj=0;jj<5;jj++) {
        comp_y[jj] = 4.0*pi/comp_y.size()*(jj+1);
        comp_data_2d[tt][jj][ii] = 0.1 * cos(comp_x[ii]+tt*dt) * sin(comp_y[jj]+tt*dt);
        comp_index_2d[tt][jj][ii] = ii + jj*100 + 10000*tt;
        for (decltype(comp_z)::size_type kk=0;kk<2;kk++) {
          comp_z[kk] = 100*(kk+1);
          comp_data_3d[tt][kk][jj][ii] = 0.1 * cos(comp_x[ii]+tt*dt) * sin(comp_y[jj]+tt*dt) + comp_z[kk];
          comp_index_3d[tt][kk][jj][ii] = ii + jj*100 + 1000*kk + 10000*tt;
        } //comp_z
      } //comp_y
    } //comp_x
  } //tt

  // Degrees of Freedom decomposition of input arrays
  std::array<Int,1> dof_x;
  std::array<Int,1> dof_y;
  std::array<Int,1> dof_z;
  std::array<Int,1> dof_test_1d;
  std::array<Int,1> dof_test_2d;
  std::array<Int,1> dof_test_3d;
  Int xstart,xstop;
  Int ystart,ystop;
  Int zstart,zstop;
  Int test_1d_start,test_1d_stop;
  Int test_2d_start,test_2d_stop;
  Int test_3d_start,test_3d_stop;
  Int myrank, numranks;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &numranks);

  get_dof(ekat::util::size(x_data), myrank, numranks, dof_x[0], xstart,xstop);
  get_dof(ekat::util::size(y_data), myrank, numranks, dof_y[0], ystart,ystop);
  get_dof(ekat::util::size(z_data), myrank, numranks, dof_z[0], zstart,zstop);
  get_dof(ekat::util::size(test_index_1d), myrank, numranks, dof_test_1d[0], test_1d_start,test_1d_stop);
  get_dof(ekat::util::size(test_index_2d), myrank, numranks, dof_test_2d[0], test_2d_start,test_2d_stop);
  get_dof(ekat::util::size(test_index_3d), myrank, numranks, dof_test_3d[0], test_3d_start,test_3d_stop);

  // Read input data and compare
  grid_read_data_array(infilename,"x",dof_x,ekat::util::data(x_data)+xstart);
  grid_read_data_array(infilename,"y",dof_y,ekat::util::data(y_data)+ystart);
  grid_read_data_array(infilename,"z",dof_z,ekat::util::data(z_data)+zstart);

  for (Int ii=xstart;ii<=xstop;ii++) {
    REQUIRE( x_data[ii] == comp_x[ii] );//2.0*pi/x_data.size()*(ii+1) );
  }
  for (Int jj=ystart;jj<=ystop;jj++) {
    REQUIRE( y_data[jj] == comp_y[jj] );//4.0*pi/y_data.size()*(jj+1) );
  }
  for (Int kk=zstart;kk<=zstop;kk++) {
    REQUIRE( z_data[kk] == comp_z[kk] );//100*(kk+1) );
  }
// TODO Loop over more timesteps
  for (int tt=0;tt<3;tt++) {
    pio_update_time(infilename,-999.0);
    grid_read_data_array(infilename,"index_1d",dof_test_1d,ekat::util::data(test_index_1d)+test_1d_start);
    grid_read_data_array(infilename,"index_2d",dof_test_2d,ekat::util::data(test_index_2d)+test_2d_start);
    grid_read_data_array(infilename,"index_3d",dof_test_3d,ekat::util::data(test_index_3d)+test_3d_start);
    grid_read_data_array(infilename,"data_1d",dof_test_1d,ekat::util::data(test_data_1d)+test_1d_start);
    grid_read_data_array(infilename,"data_2d",dof_test_2d,ekat::util::data(test_data_2d)+test_2d_start);
    grid_read_data_array(infilename,"data_3d",dof_test_3d,ekat::util::data(test_data_3d)+test_3d_start);
    for (int ii=test_1d_start;ii<test_1d_stop;ii++) {
      REQUIRE(*(ekat::util::data(test_index_1d) + ii) == *(ekat::util::data(comp_index_1d[tt]) + ii));
      REQUIRE(*(ekat::util::data(test_data_1d) + ii) == *(ekat::util::data(comp_data_1d[tt]) + ii));
    }
    for (int ii=test_2d_start;ii<test_2d_stop;ii++) {
      REQUIRE(*(ekat::util::data(test_index_2d) + ii) == *(ekat::util::data(comp_index_2d[tt]) + ii));
      REQUIRE(*(ekat::util::data(test_data_2d) + ii) == *(ekat::util::data(comp_data_2d[tt]) + ii));
    }
    for (int ii=test_3d_start;ii<test_3d_stop;ii++) {
      REQUIRE(*(ekat::util::data(test_index_3d) + ii) == *(ekat::util::data(comp_index_3d[tt]) + ii));
      REQUIRE(*(ekat::util::data(test_data_3d) + ii) == *(ekat::util::data(comp_data_3d[tt]) + ii));
    }
  } //tt

  eam_pio_finalize();
} // TEST scorpio_interface_input
/* ================================================================================================================ */
/*                                   Local functions to be used for tests:                                          */
Real f_x(const Real x, const Real t) {
  Real f;
  f = 0.1 * cos(x+t);
  return f;
}
Real f_y(const Real y, const Real t) {
  Real f;
  f = sin(y+t);
  return f;
}
Real f_z(const Real z, const Real t) {
  Real f;
  f = z;
  return f;
}
Int ind_x(const Int ii) {
  return ii;
}
Int ind_y(const Int jj) {
  return jj*100;
}
Int ind_z(const Int kk) {
  return 1000*kk;
}
Int ind_t(const Int tt) {
  return 10000*tt;
}

void get_dof(const std::size_t dimsize, const Int myrank, const Int numranks, Int dof_len, Int &istart, Int &istop) {

  dof_len = dimsize/numranks;
  Int extra_procs;
  extra_procs = dimsize % numranks;
  if (extra_procs>0) {
    dof_len = dof_len+1;
  }
  istart = myrank*dof_len;
  if (myrank == numranks-1) {
    dof_len = dimsize-istart;
  }
  istop = istart + dof_len-1;
//  std::printf("ASD MPI ranks: %d/%d : %d/%d, %d -- %d:%d\n",myrank+1,numranks,dof_len, (int) dimsize,extra_procs,istart,istop);

  return;
}
/* ================================================================================================================ */
} //namespace
