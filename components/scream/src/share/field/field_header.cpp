#include "share/field/field_header.hpp"

namespace scream
{

FieldHeader::FieldHeader (const identifier_type& id)
 : m_identifier (id)
 , m_tracking (create_tracking(id.name()))
 , m_alloc_prop ()
{
  // Nothing to be done here
}

FieldHeader::FieldHeader (const identifier_type& id,
                          std::shared_ptr<FieldHeader> parent,
                          const int idim, const int k)
 : m_identifier (id)
{
  EKAT_REQUIRE_MSG (parent!=nullptr,
      "Error! Invalid pointer for parent header.\n");
  EKAT_REQUIRE_MSG (id.get_layout_ptr()!=nullptr,
      "Error! Input field identifier has an invalid layout pointer.\n");

  m_parent = parent;

  m_tracking = create_tracking(id.name(),parent->get_tracking_ptr());

  m_alloc_prop = parent->get_alloc_properties().subview(idim,k);
  m_alloc_prop.commit(id.get_layout_ptr());

  m_tracking->register_as_children_in_parent();
}

} // namespace scream
