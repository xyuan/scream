#include "scream_scorpio_interface.hpp"
#include "ekat/ekat_scalar_traits.hpp"
#include "scream_config.h"

#include "ekat/ekat_assert.hpp"
#include "share/scream_types.hpp"

#include "gptl.h"

#include <string>

using scream::Real;
using scream::Int;
extern "C" {

// Fortran routines to be called from C++
  void register_infile_c2f(const char*&& filename);
  void set_decomp_c2f(const char*&& filename);
  void set_dof_c2f(const char*&& filename,const char*&& varname,const Int dof_len,const Int *x_dof);
  void grid_read_data_array_c2f_real(const char*&& filename, const char*&& varname, const Int dim1_length, Real *hbuf);
  void grid_read_data_array_c2f_int(const char*&& filename, const char*&& varname, const Int dim1_length, Int *hbuf);

  void grid_write_data_array_c2f_real_1d(const char*&& filename, const char*&& varname, const Int dim1_length, const Real* hbuf);
  void grid_write_data_array_c2f_real_2d(const char*&& filename, const char*&& varname, const Int dim1_length, const Int dim2_length, const Real* hbuf);
  void grid_write_data_array_c2f_real_3d(const char*&& filename, const char*&& varname, const Int dim1_length, const Int dim2_length, const Int dim3_length, const Real* hbuf);
  void grid_write_data_array_c2f_real_4d(const char*&& filename, const char*&& varname, const Int dim1_length, const Int dim2_length, const Int dim3_length, const Int dim4_length, const Real* hbuf);
  void grid_write_data_array_c2f_int_1d (const char*&& filename, const char*&& varname, const Int dim1_length, const Int* hbuf);
  void grid_write_data_array_c2f_int_2d (const char*&& filename, const char*&& varname, const Int dim1_length, const Int dim2_length, const Int* hbuf);
  void grid_write_data_array_c2f_int_3d (const char*&& filename, const char*&& varname, const Int dim1_length, const Int dim2_length, const Int dim3_length, const Int* hbuf);
  void grid_write_data_array_c2f_int_4d (const char*&& filename, const char*&& varname, const Int dim1_length, const Int dim2_length, const Int dim3_length, const Int dim4_length, const Int* hbuf);
  void eam_init_pio_subsystem_c2f(const int mpicom, const int compid, const bool local);
  void eam_pio_finalize_c2f();
  void register_outfile_c2f(const char*&& filename);
  void sync_outfile_c2f(const char*&& filename);
  void eam_pio_closefile_c2f(const char*&& filename);
  void pio_update_time_c2f(const char*&& filename,const Real time);
  void register_dimension_c2f(const char*&& filename, const char*&& shortname, const char*&& longname, const int length);
  void register_variable_c2f(const char*&& filename,const char*&& shortname, const char*&& longname, const int numdims, const char** var_dimensions, const int dtype, const char*&& pio_decomp_tag);
  void get_variable_c2f(const char*&& filename,const char*&& shortname, const char*&& longname, const int numdims, const char** var_dimensions, const int dtype, const char*&& pio_decomp_tag);
  void eam_pio_enddef_c2f(const char*&& filename);

  void count_pio_atm_file_c2f();
} // extern C

namespace scream {
namespace scorpio {
/* ----------------------------------------------------------------- */
void eam_init_pio_subsystem(const int mpicom) {
  // TODO: Right now the compid has been hardcoded to 0 and the flag
  // to create a init a subsystem in SCREAM is hardcoded to true.
  // When surface coupling is established we will need to refactor this
  // routine to pass the appropriate values depending on if we are running
  // the full model or a unit test.
  GPTLinitialize();
  eam_init_pio_subsystem_c2f(mpicom,0,true);
}
/* ----------------------------------------------------------------- */
void eam_pio_finalize() {
  eam_pio_finalize_c2f();
  GPTLfinalize();
}
/* ----------------------------------------------------------------- */
void register_outfile(const std::string& filename) {


  register_outfile_c2f(filename.c_str());
}
/* ----------------------------------------------------------------- */
void eam_pio_closefile(const std::string& filename) {

  eam_pio_closefile_c2f(filename.c_str());
}
/* ----------------------------------------------------------------- */
void register_infile(const std::string& filename) {

  register_infile_c2f(filename.c_str());
}
/* ----------------------------------------------------------------- */
void sync_outfile(const std::string& filename) {

  sync_outfile_c2f(filename.c_str());
}
/* ----------------------------------------------------------------- */
void set_decomp(const std::string& filename) {

  set_decomp_c2f(filename.c_str());
}
/* ----------------------------------------------------------------- */
void set_dof(const std::string& filename, const std::string& varname, const Int dof_len, const Int* x_dof) {

  set_dof_c2f(filename.c_str(),varname.c_str(),dof_len,x_dof);
}
/* ----------------------------------------------------------------- */
void pio_update_time(const std::string& filename, const Real time) {

  pio_update_time_c2f(filename.c_str(),time);
}
/* ----------------------------------------------------------------- */
void register_dimension(const std::string &filename, const std::string& shortname, const std::string& longname, const int length) {

  register_dimension_c2f(filename.c_str(), shortname.c_str(), longname.c_str(), length);
}
/* ----------------------------------------------------------------- */
void get_variable(const std::string &filename, const std::string& shortname, const std::string& longname, const int numdims, const std::vector<std::string>& var_dimensions, const int dtype, const std::string& pio_decomp_tag) {

  /* Convert the vector of strings that contains the variable dimensions to a char array */
  const char** var_dimensions_c = new const char*[numdims];
  for (int ii = 0;ii<numdims;++ii) 
  {
    var_dimensions_c[ii] = var_dimensions[ii].c_str();
  }
  get_variable_c2f(filename.c_str(), shortname.c_str(), longname.c_str(), numdims, var_dimensions_c, dtype, pio_decomp_tag.c_str());
  delete[] var_dimensions_c;
}
/* ----------------------------------------------------------------- */
void get_variable(const std::string &filename, const std::string& shortname, const std::string& longname, const int numdims, const char**&& var_dimensions, const int dtype, const std::string& pio_decomp_tag) {

  get_variable_c2f(filename.c_str(), shortname.c_str(), longname.c_str(), numdims, var_dimensions, dtype, pio_decomp_tag.c_str());
}
/* ----------------------------------------------------------------- */
void register_variable(const std::string &filename, const std::string& shortname, const std::string& longname, const int numdims, const std::vector<std::string>& var_dimensions, const int dtype, const std::string& pio_decomp_tag) {

  /* Convert the vector of strings that contains the variable dimensions to a char array */
  const char** var_dimensions_c = new const char*[numdims];
  for (int ii = 0;ii<numdims;++ii) 
  {
    var_dimensions_c[ii] = var_dimensions[ii].c_str();
  }
  register_variable_c2f(filename.c_str(), shortname.c_str(), longname.c_str(), numdims, var_dimensions_c, dtype, pio_decomp_tag.c_str());
  delete[] var_dimensions_c;
}
/* ----------------------------------------------------------------- */
void register_variable(const std::string &filename, const std::string& shortname, const std::string& longname, const int numdims, const char**&& var_dimensions, const int dtype, const std::string& pio_decomp_tag) {

  register_variable_c2f(filename.c_str(), shortname.c_str(), longname.c_str(), numdims, var_dimensions, dtype, pio_decomp_tag.c_str());
}
/* ----------------------------------------------------------------- */
void eam_pio_enddef(const std::string &filename) {
  eam_pio_enddef_c2f(filename.c_str());
}
/* ----------------------------------------------------------------- */
void count_pio_atm_file() {

  count_pio_atm_file_c2f();

}
/* ----------------------------------------------------------------- */
// Read and Write routines
// Inputs:
//   filename:   is the filename for the associated netCDF file.
//   varname:    is the variable name for this netCDF operation.
//   dims:       is a vector of integers with length equal to the total number of dimensions for the variable.
//               note, each dimension here is the actual physical dimension, which can be different than the allocated dimension.
//               If, for example, a variable is padded then the final value of dims would be smaller than the allocation size.
//   dim_length: is the total number of real data values.  Again, if a variable is padded, the padding should not be included in this value.
//   padding:    is the amount of padding the variable has, see documentation on packed variables in fields.hpp.  Note padding should be >= 0.
//   hbuf:       is a pointer to where the data for this variable is stored.
/* ----------------------------------------------------------------- */
// Handling the reading of input for packed arrays
void grid_read_data_array(const std::string &filename, const std::string &varname, const std::vector<int>& dims, const Int& dim_length, const Int& padding, Real *hbuf) {

  if (padding == 0) // then no xtra data, is contiguous
  {
    grid_read_data_array_c2f_real(filename.c_str(),varname.c_str(),dim_length,hbuf);
  } else {
    std::vector<Real> hbuf_new(dim_length);
    grid_read_data_array_c2f_real(filename.c_str(),varname.c_str(),dim_length,hbuf_new.data());
    // Copy the read in values back to the packed array, add padding
    add_remove_padding(dims,padding,hbuf_new.data(),hbuf,true);
  }
};
/* ----------------------------------------------------------------- */
void grid_read_data_array(const std::string &filename, const std::string &varname, const Int& dim_length, Real *hbuf) {

  grid_read_data_array_c2f_real(filename.c_str(),varname.c_str(),dim_length,hbuf);

};
/* ----------------------------------------------------------------- */
// Handling the writing of output for packed arrays
void grid_write_data_array(const std::string &filename, const std::string &varname, const std::vector<int>& dims, const Int& dim_length, const Int& padding, const Real* hbuf)
{
  if (padding == 0) // then no xtra data, is contiguous
  {
    grid_write_data_array_c2f_real_1d(filename.c_str(),varname.c_str(),dim_length,hbuf);
  } else {
    std::vector<Real> hbuf_new(dim_length);
    // Packed along final dimension
    add_remove_padding(dims,padding,hbuf,hbuf_new.data(),false);
    grid_write_data_array_c2f_real_1d(filename.c_str(),varname.c_str(),dim_length,hbuf_new.data());
  }
}
/* ----------------------------------------------------------------- */
void grid_write_data_array(const std::string &filename, const std::string &varname, const Int& dim_length, const Real* hbuf) {

  grid_write_data_array_c2f_real_1d(filename.c_str(),varname.c_str(),dim_length,hbuf);

};
/* ----------------------------------------------------------------- */
void add_remove_padding(const std::vector<int>& dims, const int padding,
                        const Real* data_in, Real* data_out, const bool add_padding)
{
  using KT = KokkosTypes<HostDevice>;

  // Split dims in [dim_1,...,dim_{N-1}] and dim_N, compute product P of first N-1 dims,
  // and reshape pointers as 2d views of dims (P,dim_N) and (P,dim_N+padding)
  int fast_dim = dims.back();
  int size = 1;
  for (auto d : dims) {
    size *= d;
  }
  int lumped_slow_dims = size / fast_dim;

  const int dim2_in  = add_padding ? fast_dim : fast_dim + padding;
  const int dim2_out = add_padding ? fast_dim + padding : fast_dim;
  ekat::Unmanaged<KT::view_2d<const Real>> from(data_in, lumped_slow_dims,dim2_in);
  ekat::Unmanaged<KT::view_2d<Real>>       to  (data_out,lumped_slow_dims,dim2_out);
  for (int i=0; i<lumped_slow_dims; ++i) {
    for (int j=0; j<fast_dim; ++j) {
      to(i,j) = from(i,j);
    }
    if (add_padding) {
      for (int j=fast_dim; j<dim2_out; ++j) {
        to(i,j) = ekat::ScalarTraits<Real>::invalid();
      }
    }
  }
}
/* ----------------------------------------------------------------- */

} // namespace scorpio
} // namespace scream
