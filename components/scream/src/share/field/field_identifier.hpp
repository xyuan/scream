#ifndef SCREAM_FIELD_IDENTIFIER_HPP
#define SCREAM_FIELD_IDENTIFIER_HPP

#include "share/field/field_layout.hpp"
#include "ekat/util/ekat_string_utils.hpp"
#include "ekat/util/ekat_units.hpp"
#include "ekat/std_meta/ekat_std_enable_shared_from_this.hpp"

#include <vector>

namespace scream
{

/*
 *  A small class to hold basic info about a field
 *
 *  This tiny class is only used to uniquely identify a field, using
 *  the name, and its layout (which containss the rank, the tags,
 *  and possibly the dimensions).
 *  There is no additional meta data about this field.
 */

class FieldIdentifier {
public:
  using layout_type     = FieldLayout;
  using layout_ptr_type = std::shared_ptr<const layout_type>;
  using ci_string       = ekat::CaseInsensitiveString;
  using units_type      = ekat::units::Units;

  // Constructor(s)
  FieldIdentifier () = delete;
  FieldIdentifier (const FieldIdentifier&) = default;
  FieldIdentifier (const std::string& name,
                   const layout_type& layout,
                   const units_type& units,
                   const std::string& grid_name);

  // Delete assignment, to prevent overwriting identifiers sneakyly
  FieldIdentifier& operator= (const FieldIdentifier&) = delete;

  // ----- Getters ----- //

  // Name and layout informations
  const std::string&      name           () const { return  m_name;      }
  const layout_type&      get_layout     () const { return *m_layout;    }
  const layout_ptr_type&  get_layout_ptr () const { return  m_layout;    }
  const units_type&       get_units      () const { return  m_units;     }
  const std::string&      get_grid_name  () const { return  m_grid_name; }

  // The identifier string is a conveniet way to display the information of
  // the identifier, so that it can be easily read.
  std::string get_id_string () const;

  // ----- Setters ----- //

  // Note: as soon as the layout is set, it cannot be changed.
  void set_layout (const layout_type& layout);
  void set_layout (const layout_ptr_type& layout);

  // We reimplement the equality operator for identifiers comparison (needed for some std container)
  friend bool operator== (const FieldIdentifier&, const FieldIdentifier&);
  friend bool operator<  (const FieldIdentifier&, const FieldIdentifier&);

protected:

  // The field name
  ci_string       m_name;

  // The layout of the field
  layout_ptr_type m_layout;

  // The units of this field
  units_type      m_units;

  // The name of the grid on which the field is defined
  ci_string       m_grid_name;
};

bool operator== (const FieldIdentifier& fid1, const FieldIdentifier& fid2);
inline bool operator!= (const FieldIdentifier& fid1, const FieldIdentifier& fid2) { return !(fid1==fid2); }

} // namespace scream

#endif // SCREAM_FIELD_IDENTIFIER_HPP
