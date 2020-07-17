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
void get_dof(const std::size_t dimsize, const Int myrank, const Int numranks, Int &dof_len, Int &istart, Int &istop);

TEST_CASE("scorpio_interface_output", "") {
/* A test to check that both output and input are behaving properly in the scream interface.
 * The first step creates a new output file, generates data for writing output, and commits that data to the new output file.
 * The second step opens the recently created output file, reads in all the available data, 
 * and then compares the read in values with the expected values that should have been written out.
 */

  /* load namespaces needed for I/O and for data management */
  using namespace scream;
  using namespace scream::scorpio;
  using ekat::util::data;

  int nerr = 0; // Set a record of the number of errors encountered.
  /* Create the set of SCORPIO output files and their respective dimensions and variables. */
  int compid=0;  // For CIME based builds this will be the integer ID assigned to the atm by the component coupler.  For testing we simply set to 0
  Int myrank, numranks;
  MPI_Fint fcomm = MPI_Comm_c2f(MPI_COMM_WORLD);  // MPI communicator group used for I/O.  In our simple test we use MPI_COMM_WORLD, however a sub
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);  // Store rank and total number of ranks for determining which chunk of global array this rank is responsible for reading.
  MPI_Comm_size(MPI_COMM_WORLD, &numranks);
  eam_init_pio_subsystem(fcomm,compid,true);   // Gather the initial PIO subsystem data creater by component coupler
  /* Tell the scorpio interface module that we have a new output file to write to */
  std::string outfilename = "scorpio_output_test.nc";  // For simplicity create a variable to store the name of the output file.  This will be used with most calls to make sure the every function is applying the action to the proper file (if multiple files are open)
  register_outfile(outfilename);                       // Register this output file with the scorpio interface module
  /* Set up the 4 spatial dimensions for this test and register them with the new output file */
  int xlen=10, ylen=5, zlen=2;
  register_dimension(outfilename,"x","horizontal distance",xlen);
  register_dimension(outfilename,"y","vertical distance",ylen);
  register_dimension(outfilename,"z","height",zlen);
  register_dimension(outfilename,"time","time",0);  // Note that time has an unknown length, setting the "length" to 0 tells the interface to set this dimension as having an unlimited length, thus allowing us to write as many timesnaps to file as we desire.
  /* Set up the list variables we wish to include in the output file
   * Here we note that we are creating character arrays which use the named dimensions above.
   * This is used when we register the variables to tell scorpio what dimensions each variable uses.
   */
  const char* vec_time[] = {"time"};
  const char* vec_x[]    = {"x"};
  const char* vec_y[]    = {"y"};
  const char* vec_z[]    = {"z"};
  const char* vec_xt[]   = {"x","time"};
  const char* vec_xyt[]  = {"x","y","time"};
  const char* vec_xyzt[] = {"x","y","z","time"};
  const char* vec_xy[]   = {"x","y"}; 
  const char* vec_xyz[]  = {"x","y","z"};
 
  register_variable(outfilename,"time","time",1,vec_time,  PIO_REAL,"t");
  register_variable(outfilename,"x","x-direction",1,vec_x, PIO_REAL,"x-real");
  register_variable(outfilename,"y","y-direction",1,vec_y, PIO_REAL,"y-real");
  register_variable(outfilename,"z","z-direction",1,vec_z, PIO_REAL,"z-real");
  register_variable(outfilename,"data_1d","test value for 1d field",2,vec_xt,   PIO_REAL,"xt-real");
  register_variable(outfilename,"data_2d","test value for 2d field",3,vec_xyt,  PIO_REAL,"xyt-real");
  register_variable(outfilename,"data_3d","test value for 3d field",4,vec_xyzt, PIO_REAL,"xyzt-real");
  register_variable(outfilename,"index_1d","test value for 1d field",2,vec_xt,   PIO_INT,"xt-int");
  register_variable(outfilename,"index_2d","test value for 2d field",3,vec_xyt,  PIO_INT,"xyt-int");
  register_variable(outfilename,"index_3d","test value for 3d field",4,vec_xyzt, PIO_INT,"xyzt-int");
  /* When we are finished defining the set of dimensions and variables for the new output file
   * we have to officially "end" the definition phase in scorpio.  After this step we can no longer
   * add dimensions or variables to the file.
   */
  eam_pio_enddef(outfilename);
  /* Construst the data to be written as output.  Note that here we take advantage of the ekat
   * md_array utility for multi-dimension arrays.  It would also be acceptable to use the std::array,
   * or other construct.  In the end, what is passed to scorpio is only the pointer to the beginning
   * of the array of data.
   */
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
  /* Set the x, y and z dimension values */
  for (decltype(x_data)::size_type ii=0;ii<x_data.size();ii++) {
    x_data[ii] = 2.0*pi/x_data.size()*(ii+1);
  }
  for (int jj=0;jj<5;jj++) {
    y_data[jj] = 4.0*pi/y_data.size()*(jj+1);
  }
  for (int kk=0;kk<2;kk++) {
    z_data[kk] = 100*(kk+1);
  }
  /* Write these values to the file.  Note, "x", "y" and "z" are both
   * dimensions and variables.  Here we are writing to the "variable"
   * definition.  Dimensions in scorpio don't have vectors associated
   * with them, they are labels with lengths used to define variables.
   * So here we note, that you will most likely want to define all
   * dimensions as variables as well.
   */
  grid_write_data_array(outfilename,"x",xdim,ekat::util::data(x_data));
  grid_write_data_array(outfilename,"y",ydim,ekat::util::data(y_data));
  grid_write_data_array(outfilename,"z",zdim,ekat::util::data(z_data));
//  sync_outfile(outfilename); 
  /* To test multiple timelevels, generate unique data for multiple timesnaps, then write the data to file. */
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
  /* Now close the output file and reopen it to check that the output is correct. */
  eam_pio_closefile(outfilename);
  register_infile(outfilename);
  //
  get_variable(outfilename,"time","time",1,vec_time, PIO_REAL,"t");
  get_variable(outfilename,"x","x-direction",1,vec_x, PIO_REAL,"x-real");
  get_variable(outfilename,"y","y-direction",1,vec_y, PIO_REAL,"y-real");
  get_variable(outfilename,"z","z-direction",1,vec_z, PIO_REAL,"z-real");
  get_variable(outfilename,"data_1d","test value for 1d field",1,vec_x, PIO_REAL,"x-real");
  get_variable(outfilename,"data_2d","test value for 2d field",2,vec_xy, PIO_REAL,"xy-real");
  get_variable(outfilename,"data_3d","test value for 3d field",3,vec_xyz, PIO_REAL,"xyz-real");
  get_variable(outfilename,"index_1d","test value for 1d field",1,vec_x, PIO_INT,"x-int");
  get_variable(outfilename,"index_2d","test value for 2d field",2,vec_xy, PIO_INT,"xy-int");
  get_variable(outfilename,"index_3d","test value for 3d field",3,vec_xyz, PIO_INT,"xyz-int");
  set_decomp(outfilename);
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
  // Setup PIO configuation for reading data
  get_dof(ekat::util::size(x_data), myrank, numranks, dof_x[0], xstart,xstop);
  std::vector<Int> x_dof(dof_x[0]);
  for (int ii=xstart,cnt=0;ii<=xstop;++ii,++cnt) {
    x_dof[cnt] = ii;
  }
  set_dof(outfilename,"x",dof_x[0],x_dof.data());

  get_dof(ekat::util::size(y_data), myrank, numranks, dof_y[0], ystart,ystop);
  std::vector<Int> y_dof(dof_y[0]);
  for (int ii=ystart,cnt=0;ii<=ystop;++ii,++cnt) {
    y_dof[cnt] = ii;
  }
  set_dof(outfilename,"y",dof_y[0],y_dof.data());

  get_dof(ekat::util::size(z_data), myrank, numranks, dof_z[0], zstart,zstop);
  std::vector<Int> z_dof(dof_z[0]);
  for (int ii=zstart,cnt=0;ii<=zstop;++ii,++cnt) {
    z_dof[cnt] = ii;
  }
  set_dof(outfilename,"z",dof_z[0],z_dof.data());

  get_dof(dimlen_1d[0], myrank, numranks, dof_test_1d[0], test_1d_start,test_1d_stop);
  std::vector<Int> test1d_dof(dof_test_1d[0]);
  for (int ii=test_1d_start,cnt=0;ii<=test_1d_stop;++ii,++cnt) {
    test1d_dof[cnt] = ii;
  }
  set_dof(outfilename,"index_1d",dof_test_1d[0],test1d_dof.data());
  set_dof(outfilename,"data_1d",dof_test_1d[0],test1d_dof.data());

  get_dof(dimlen_2d[0]*dimlen_2d[1], myrank, numranks, dof_test_2d[0], test_2d_start,test_2d_stop);
  std::vector<Int> test2d_dof(dof_test_2d[0]);
  for (int ii=test_2d_start,cnt=0;ii<=test_2d_stop;++ii,++cnt) {
    test2d_dof[cnt] = ii;
  }
  set_dof(outfilename,"index_2d",dof_test_2d[0],test2d_dof.data());
  set_dof(outfilename,"data_2d",dof_test_2d[0],test2d_dof.data());

  get_dof(dimlen_3d[0]*dimlen_3d[1]*dimlen_3d[2], myrank, numranks, dof_test_3d[0], test_3d_start,test_3d_stop);
  std::vector<Int> test3d_dof(dof_test_3d[0]);
  for (int ii=test_3d_start,cnt=0;ii<=test_3d_stop;++ii,++cnt) {
    test3d_dof[cnt] = ii;
  }
  set_dof(outfilename,"index_3d",dof_test_3d[0],test3d_dof.data());
  set_dof(outfilename,"data_3d",dof_test_3d[0],test3d_dof.data());
  // Read input data and compare
  grid_read_data_array(outfilename,"x",dof_x[0],ekat::util::data(x_data)+xstart);
  grid_read_data_array(outfilename,"y",dof_y[0],ekat::util::data(y_data)+ystart);
  grid_read_data_array(outfilename,"z",dof_z[0],ekat::util::data(z_data)+zstart);

  for (Int ii=xstart;ii<=xstop;ii++) {
    if (x_data[ii] != 2.0*pi/x_data.size()*(ii+1)) { ++nerr;}
  }
  REQUIRE(nerr==0);
  nerr = 0;
  for (Int jj=ystart;jj<=ystop;jj++) {
    if (y_data[jj] != 4.0*pi/y_data.size()*(jj+1)) { ++nerr;}
  }
  REQUIRE(nerr==0);
  nerr = 0;
  for (Int kk=zstart;kk<=zstop;kk++) {
    if (z_data[kk] != 100*(kk+1)) { ++nerr;}
  }
  REQUIRE(nerr==0);
  nerr = 0;
  for (int tt=0;tt<3;tt++) {
    pio_update_time(outfilename,-999.0);
    grid_read_data_array(outfilename,"index_1d",dof_test_1d[0],ekat::util::data(test_index_1d)+test_1d_start);
    grid_read_data_array(outfilename,"index_2d",dof_test_2d[0],ekat::util::data(test_index_2d)+test_2d_start);
    grid_read_data_array(outfilename,"index_3d",dof_test_3d[0],ekat::util::data(test_index_3d)+test_3d_start);
    grid_read_data_array(outfilename,"data_1d", dof_test_1d[0],ekat::util::data(test_data_1d) +test_1d_start);
    grid_read_data_array(outfilename,"data_2d", dof_test_2d[0],ekat::util::data(test_data_2d) +test_2d_start);
    grid_read_data_array(outfilename,"data_3d", dof_test_3d[0],ekat::util::data(test_data_3d) +test_3d_start);
    nerr = 0;
    for (int ii=0,ind=test_1d_start;ii<x_data.size();ii++,ind++) {
      if (*(ekat::util::data(test_data_1d) + ind)  != f_x(x_data[ii],tt*dt)) {++nerr;}
      if (*(ekat::util::data(test_index_1d) + ind) != ind_x(ii) + ind_t(tt))  {++nerr;}
    }
    REQUIRE(nerr==0);
    nerr = 0;
    for (int jj=0,ind=test_2d_start;jj<y_data.size();jj++) {
      for (int ii=0;ii<x_data.size();ii++,ind++) {
        if (*(ekat::util::data(test_data_2d) + ind)  != f_x(x_data[ii],tt*dt)*f_y(y_data[jj],tt*dt)) {++nerr;}
        if (*(ekat::util::data(test_index_2d) + ind) != ind_y(jj) + ind_x(ii) + ind_t(tt))  {++nerr;}
      }
    }
    REQUIRE(nerr==0);
    nerr = 0;
    for (int kk=0,ind=test_3d_start;kk<z_data.size();kk++) {
      for (int jj=0;jj<y_data.size();jj++) {
        for (int ii=0;ii<x_data.size();ii++,ind++) {
          if (*(ekat::util::data(test_data_3d) + ind)  != f_x(x_data[ii],tt*dt)*f_y(y_data[jj],tt*dt) + f_z(z_data[kk],tt*dt)) {++nerr;}
          if (*(ekat::util::data(test_index_3d) + ind) != ind_z(kk) + ind_y(jj) + ind_x(ii) + ind_t(tt))  {++nerr;}
        }
      }
    }
    REQUIRE(nerr==0);
  } //tt

  eam_pio_finalize();
} // TEST scorpio_interface_output
/* ================================================================================================================ */
TEST_CASE("scorpio_interface_input", "") {

  using namespace scream;
  using namespace scream::scorpio;
  using ekat::util::data;

  Real pi = 2*acos(0.0);
  Real dt = 1.0;
  int nerr = 0;
  int compid=0;
  Int myrank, numranks;
  MPI_Fint fcomm = MPI_Comm_c2f(MPI_COMM_WORLD);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &numranks);
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
 
  get_variable(infilename,"time","time",1,vec_time, PIO_REAL,"t");
  get_variable(infilename,"x","x-direction",1,vec_x, PIO_REAL,"x-real");
  get_variable(infilename,"y","y-direction",1,vec_y, PIO_REAL,"y-real");
  get_variable(infilename,"z","z-direction",1,vec_z, PIO_REAL,"z-real");
  get_variable(infilename,"data_1d","test value for 1d field",1,vec_x, PIO_REAL,"x-real");
  get_variable(infilename,"data_2d","test value for 2d field",2,vec_xy, PIO_REAL,"xy-real");
  get_variable(infilename,"data_3d","test value for 3d field",3,vec_xyz, PIO_REAL,"xyz-real");
  get_variable(infilename,"index_1d","test value for 1d field",1,vec_x, PIO_INT,"x-int");
  get_variable(infilename,"index_2d","test value for 2d field",2,vec_xy, PIO_INT,"xy-int");
  get_variable(infilename,"index_3d","test value for 3d field",3,vec_xyz, PIO_INT,"xyz-int");
  set_decomp(infilename);

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
  ekat::util::md_array<Real,3,10>       comp_data_1d;
  ekat::util::md_array<Real,3, 5,10>    comp_data_2d;
  ekat::util::md_array<Real,3, 2, 5,10> comp_data_3d;
  ekat::util::md_array<Int, 3,10>       comp_index_1d;
  ekat::util::md_array<Int, 3, 5,10>    comp_index_2d;
  ekat::util::md_array<Int, 3, 2, 5,10> comp_index_3d;
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

  get_dof(ekat::util::size(x_data), myrank, numranks, dof_x[0], xstart,xstop);
  std::vector<Int> x_dof(dof_x[0]);
  for (int ii=xstart,cnt=0;ii<=xstop;++ii,++cnt) {
    x_dof[cnt] = ii;
  }
  set_dof(infilename,"x",dof_x[0],x_dof.data());

  get_dof(ekat::util::size(y_data), myrank, numranks, dof_y[0], ystart,ystop);
  std::vector<Int> y_dof(dof_y[0]);
  for (int ii=ystart,cnt=0;ii<=ystop;++ii,++cnt) {
    y_dof[cnt] = ii;
  }
  set_dof(infilename,"y",dof_y[0],y_dof.data());

  get_dof(ekat::util::size(z_data), myrank, numranks, dof_z[0], zstart,zstop);
  std::vector<Int> z_dof(dof_z[0]);
  for (int ii=zstart,cnt=0;ii<=zstop;++ii,++cnt) {
    z_dof[cnt] = ii;
  }
  set_dof(infilename,"z",dof_z[0],z_dof.data());

  get_dof(dimlen_1d[0], myrank, numranks, dof_test_1d[0], test_1d_start,test_1d_stop);
  std::vector<Int> test1d_dof(dof_test_1d[0]);
  for (int ii=test_1d_start,cnt=0;ii<=test_1d_stop;++ii,++cnt) {
    test1d_dof[cnt] = ii;
  }
  set_dof(infilename,"index_1d",dof_test_1d[0],test1d_dof.data());
  set_dof(infilename,"data_1d",dof_test_1d[0],test1d_dof.data());

  get_dof(dimlen_2d[0]*dimlen_2d[1], myrank, numranks, dof_test_2d[0], test_2d_start,test_2d_stop);
  std::vector<Int> test2d_dof(dof_test_2d[0]);
  for (int ii=test_2d_start,cnt=0;ii<=test_2d_stop;++ii,++cnt) {
    test2d_dof[cnt] = ii;
  }
  set_dof(infilename,"index_2d",dof_test_2d[0],test2d_dof.data());
  set_dof(infilename,"data_2d",dof_test_2d[0],test2d_dof.data());

  get_dof(dimlen_3d[0]*dimlen_3d[1]*dimlen_3d[2], myrank, numranks, dof_test_3d[0], test_3d_start,test_3d_stop);
  std::vector<Int> test3d_dof(dof_test_3d[0]);
  for (int ii=test_3d_start,cnt=0;ii<=test_3d_stop;++ii,++cnt) {
    test3d_dof[cnt] = ii;
  }
  set_dof(infilename,"index_3d",dof_test_3d[0],test3d_dof.data());
  set_dof(infilename,"data_3d",dof_test_3d[0],test3d_dof.data());

  // Read input data and compare
  grid_read_data_array(infilename,"x",dof_x[0],ekat::util::data(x_data)+xstart);
  grid_read_data_array(infilename,"y",dof_y[0],ekat::util::data(y_data)+ystart);
  grid_read_data_array(infilename,"z",dof_z[0],ekat::util::data(z_data)+zstart);

  for (Int ii=xstart;ii<=xstop;ii++) {
    if (x_data[ii] != comp_x[ii]) { ++nerr;}
  }
  REQUIRE(nerr==0);
  nerr = 0;
  for (Int jj=ystart;jj<=ystop;jj++) {
    if (y_data[jj] != comp_y[jj]) { ++nerr;}
  }
  REQUIRE(nerr==0);
  nerr = 0;
  for (Int kk=zstart;kk<=zstop;kk++) {
    if (z_data[kk] != comp_z[kk]) { ++nerr;}
  }
  REQUIRE(nerr==0);
  nerr = 0;
  for (int tt=0;tt<3;tt++) {
    pio_update_time(infilename,-999.0);
    grid_read_data_array(infilename,"index_1d",dof_test_1d[0],ekat::util::data(test_index_1d)+test_1d_start);
    grid_read_data_array(infilename,"index_2d",dof_test_2d[0],ekat::util::data(test_index_2d)+test_2d_start);
    grid_read_data_array(infilename,"index_3d",dof_test_3d[0],ekat::util::data(test_index_3d)+test_3d_start);
    grid_read_data_array(infilename,"data_1d", dof_test_1d[0],ekat::util::data(test_data_1d) +test_1d_start);
    grid_read_data_array(infilename,"data_2d", dof_test_2d[0],ekat::util::data(test_data_2d) +test_2d_start);
    grid_read_data_array(infilename,"data_3d", dof_test_3d[0],ekat::util::data(test_data_3d) +test_3d_start);
    for (int ii=test_1d_start;ii<test_1d_stop;ii++) {
      if (*(ekat::util::data(test_index_1d) + ii) != *(ekat::util::data(comp_index_1d[tt]) + ii)) {++nerr;}
      if (*(ekat::util::data(test_data_1d) + ii)  != *(ekat::util::data(comp_data_1d[tt]) + ii))  {++nerr;}
    }
    REQUIRE(nerr==0);
    nerr = 0;
    for (int ii=test_2d_start;ii<test_2d_stop;ii++) {
      if (*(ekat::util::data(test_index_2d) + ii) != *(ekat::util::data(comp_index_2d[tt]) + ii)) {++nerr;}
      if (*(ekat::util::data(test_data_2d) + ii)  != *(ekat::util::data(comp_data_2d[tt]) + ii))  {++nerr;}
    }
    REQUIRE(nerr==0);
    nerr = 0;
    for (int ii=test_3d_start;ii<test_3d_stop;ii++) {
      if (*(ekat::util::data(test_index_3d) + ii) != *(ekat::util::data(comp_index_3d[tt]) + ii)) {++nerr;}
      if (*(ekat::util::data(test_data_3d) + ii)  != *(ekat::util::data(comp_data_3d[tt]) + ii))  {++nerr;}
    }
    REQUIRE(nerr==0);
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

void get_dof(const std::size_t dimsize, const Int myrank, const Int numranks, Int &dof_len, Int &istart, Int &istop) {

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

  return;
}
/* ================================================================================================================ */
} //namespace
