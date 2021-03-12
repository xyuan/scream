#include "share/field/field_identifier.hpp"

namespace scream
{

FieldIdentifier::
FieldIdentifier (const std::string& name,
                 const layout_type& layout,
                 const units_type& units,
                 const std::string& grid_name)
 : FieldIdentifier(name,"",layout,units,grid_name)
{
  // Nothing to do here
}

FieldIdentifier::
FieldIdentifier (const std::string& name,
                 const std::string& long_name,
                 const layout_type& layout,
                 const units_type&  units,
                 const std::string& grid_name)
{
  m_data = std::make_shared<Data>(name,long_name,layout,units,grid_name);
}

void FieldIdentifier::set_layout (const layout_type& layout) {
  EKAT_REQUIRE_MSG (!m_data->layout.are_dimensions_set(),
      "Error! You cannot reset the layout once it's set.\n");
  EKAT_REQUIRE_MSG (layout.are_dimensions_set(),
      "Error! Input layout must have dimensions set.\n");

  m_data->layout = layout;
}

void FieldIdentifier::set_long_name (const std::string& long_name) {
  EKAT_REQUIRE_MSG (m_data->long_name=="",
      "Error! You cannot change the long name once it has been set.\n");

  m_data->long_name = long_name;
}

std::string FieldIdentifier::get_id_string () const {
  // Create a verbose identifier string.
  std::string id = m_data->name + "[" + m_data->grid_name + "]";
  id += "<" + e2str(m_data->layout.tags()[0]);
  for (int dim=1; dim<m_data->layout.rank(); ++dim) {
    id += "," + e2str(m_data->layout.tags()[dim]);
  }
  id += ">(" + std::to_string(m_data->layout.dims()[0]);
  for (int dim=1; dim<m_data->layout.rank(); ++dim) {
    id += "," + std::to_string(m_data->layout.dims()[dim]);
  }
  id += ") [" + m_data->units.get_string() + "]";
  return id;
}

// Free functions for identifiers comparison
bool operator== (const FieldIdentifier& fid1, const FieldIdentifier& fid2) {
  // Simply compare the identifiers
  return (fid1.get_id_string()==fid2.get_id_string());
}

bool operator< (const FieldIdentifier& fid1, const FieldIdentifier& fid2) {
  // Simply compare the identifiers
  return (fid1.get_id_string()<fid2.get_id_string());
}

} // namespace scream
