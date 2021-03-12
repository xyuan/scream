#ifndef SCREAM_FIELD_IDENTIFIER_HPP
#define SCREAM_FIELD_IDENTIFIER_HPP

#include "share/field/field_layout.hpp"
#include "ekat/util/ekat_string_utils.hpp"
#include "ekat/util/ekat_units.hpp"

#include <vector>

namespace scream
{

/*
 *  A small class to hold basic info about a field
 *
 *  The FieldIdentifier (FID) is used to uniquely identify a field, using
 *  the name, its layout (tags and dimensions), the name of the grid
 *  where the field is defined, and its units.
 *  The FID stores 2 names: a short one (m_name) and a long one (m_long_name).
 *  The latter is built according to CF names guidelines:
 *    http://cfconventions.org/Data/cf-standard-names/docs/guidelines.html
 *  The user only needs to pass a name at construction time, but
 *  can also provide the long CF-compliant one (if known).
 *  When field are registered in the FieldRepo (FR), we check that both
 *  short and long names are unique, meaning that
 *   - the same short name is not used for two different long names;
 *   - the same long name is not used for two different short names.
 *  The FR will also take care of ensuring that a valid long name
 *  is stored.
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

  FieldIdentifier (const std::string& name,
                   const std::string& long_name,
                   const layout_type& layout,
                   const units_type& units,
                   const std::string& grid_name);

  // Delete assignment, to prevent overwriting identifiers sneakyly
  FieldIdentifier& operator= (const FieldIdentifier&) = delete;

  // ----- Getters ----- //

  // Name and layout informations
  const std::string&      name           () const { return  m_data->name;      }
  const std::string&      long_name      () const { return  m_data->long_name; }
  const layout_type&      get_layout     () const { return *m_data->layout;    }
  const layout_ptr_type&  get_layout_ptr () const { return  m_data->layout;    }
  const units_type&       get_units      () const { return  m_data->units;     }
  const std::string&      get_grid_name  () const { return  m_data->grid_name; }

  // The identifier string is a conveniet way to display the information of
  // the identifier, so that it can be easily read.
  std::string get_id_string () const;

  // ----- Setters ----- //

  // Note: as soon as the layout is set, it cannot be changed.
  void set_layout (const layout_type& layout);
  void set_layout (const layout_ptr_type& layout);

  // Note: as soon as the long name is set, it cannot be changed
  void set_long_name (const std::string& long_name);

  // We reimplement the equality operator for identifiers comparison (needed for some std container)
  friend bool operator== (const FieldIdentifier&, const FieldIdentifier&);
  friend bool operator<  (const FieldIdentifier&, const FieldIdentifier&);

protected:

  struct Data {
    Data (const std::string& _name,
          const std::string& _long_name,
          const layout_type& _layout,
          const units_type&  _units,
          const std::string& _grid_name)
      : name(_name)
      , long_name(_long_name)
      , layout(std::make_shared<const layout_type>(_layout))
      , units(_units)
      , grid_name(_grid_name)
    {
      // Nothing to do here
    }

    // The field name
    ci_string       name;

    // An long_name, compliant with cf naming conventions
    ci_string       long_name;

    // The layout of the field
    layout_ptr_type layout;

    // The units of this field
    units_type      units;

    // The name of the grid on which the field is defined
    ci_string       grid_name;
  };

  std::shared_ptr<Data> m_data;
};

bool operator== (const FieldIdentifier& fid1, const FieldIdentifier& fid2);
inline bool operator!= (const FieldIdentifier& fid1, const FieldIdentifier& fid2) { return !(fid1==fid2); }

} // namespace scream

#endif // SCREAM_FIELD_IDENTIFIER_HPP
