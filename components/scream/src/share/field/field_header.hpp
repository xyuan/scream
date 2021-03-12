#ifndef SCREAM_FIELD_HEADER_HPP
#define SCREAM_FIELD_HEADER_HPP

#include "share/field/field_identifier.hpp"
#include "share/field/field_tracking.hpp"
#include "share/field/field_alloc_prop.hpp"
#include "share/scream_types.hpp"

#include "share/util/scream_time_stamp.hpp"
#include "ekat/std_meta/ekat_std_any.hpp"

#include <vector>
#include <map>
#include <memory>   // For std::shared_ptr and std::weak_ptr

namespace scream
{

class AtmosphereProcess;

/*
 * A small class to contain meta-data about a field
 *
 * The FieldHeader class is itself a container of other
 * more speicific classes, such as FieldIdentifier
 * (which contains information used to uniquely identify
 * the field) or FieldTracking (which contains info used
 * to track access to the field).
 */

class FieldHeader {
public:

  using identifier_type = FieldIdentifier;
  using tracking_type   = FieldTracking;

  // Constructor(s)
  FieldHeader (const FieldHeader&) = default;
  explicit FieldHeader (const identifier_type& id);
  FieldHeader (const identifier_type& id,
               std::shared_ptr<FieldHeader> parent,
               const int idim, const int k);

  // Assignment deleted, to prevent sneaky overwrites.
  FieldHeader& operator= (const FieldHeader&) = delete;

  // ----- Getters ----- //

  // Get the basic information from the identifier
  const identifier_type& get_identifier () const { return m_identifier; }

  // Get the tracking
  const tracking_type& get_tracking () const { return *m_tracking; }
        tracking_type& get_tracking ()       { return *m_tracking; }

  // Get the allocation properties
  const FieldAllocProp& get_alloc_properties () const { return m_alloc_prop; }
        FieldAllocProp& get_alloc_properties ()       { return m_alloc_prop; }

  // Get parent (if any)
  std::weak_ptr<FieldHeader> get_parent () const { return m_parent; }

protected:

  // Static information about the field: name, rank, tags
  identifier_type                 m_identifier;

  // Tracking of the field
  std::shared_ptr<tracking_type>  m_tracking;

  // Allocation properties
  FieldAllocProp                  m_alloc_prop;

  // If this field is a sub-view of another field, we keep a pointer to the parent
  std::weak_ptr<FieldHeader>      m_parent;

};

} // namespace scream

#endif // SCREAM_FIELD_HEADER_HPP
