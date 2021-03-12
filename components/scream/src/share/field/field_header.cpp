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
  EKAT_REQUIRE_MSG (id.get_layout().are_dimensions_set(),
      "Error! Input field identifier has an incomplete layout.\n");

  m_parent = parent;

  m_tracking = create_tracking(id.name(),parent->m_tracking);

  m_alloc_prop = parent->get_alloc_properties().subview(idim,k);
  m_alloc_prop.commit(id.get_layout());

  m_tracking->register_as_children_in_parent();
}

} // namespace scream
