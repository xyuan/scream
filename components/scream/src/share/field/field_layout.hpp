#ifndef SCREAM_FIELD_LAYOUT_HPP
#define SCREAM_FIELD_LAYOUT_HPP

#include "field_tag.hpp"
#include <ekat/std_meta/ekat_std_utils.hpp>
#include <ekat/ekat_assert.hpp>

#include <string>
#include <vector>

namespace scream
{

/*
 *  A small class to hold basic info about a field layout
 *
 */

class FieldLayout {
public:

  // Constructor(s)
  FieldLayout ();
  FieldLayout (const FieldLayout&) = default;

  FieldLayout (const std::vector<FieldTag>& tags);

  FieldLayout (const std::vector<FieldTag>& tags,
               const std::vector<int>& dims);

  // Assignment (defaulted)
  FieldLayout& operator= (const FieldLayout&) = default;

  // ----- Getters ----- //

  // Name and layout informations
  const std::vector<FieldTag>& tags () const { return m_data->tags; }
  FieldTag tag  (const int idim) const;
  bool has_tag (const FieldTag t) const { return ekat::contains(m_data->tags,t); }

  // The rank is the number of tags associated to this field.
  int rank () const  { return m_data->tags.size(); }

  int dim (const FieldTag tag) const;
  int dim (const int idim) const;
  const std::vector<int>& dims () const { return m_data->dims; }

  int size ()               const;

  bool are_dimensions_set () const;

  // ----- Setters ----- //

  // Note: as soon as dimensions are set, they cannot be changed.
  void set_dimensions (const std::vector<int>& dims);

protected:

  struct Data {
    Data () = default;

    std::vector<FieldTag> tags;
    std::vector<int>      dims;
  };

  std::shared_ptr<Data> m_data;
};

bool operator== (const FieldLayout& fl1, const FieldLayout& fl2);

// ========================== IMPLEMENTATION ======================= //

inline int FieldLayout::dim (const FieldTag t) const {
  auto it = ekat::find(m_data->tags,t);

  // Check if found
  EKAT_REQUIRE_MSG(it!=m_data->tags.end(), "Error! Tag '" + e2str(t) + "' not found.\n");

  // Check only one tag (no ambiguity)
  EKAT_REQUIRE_MSG(ekat::count(m_data->tags,t)==1,
                     "Error! Tag '" + e2str(t) + "' appears multiple times.\n"
                     "       You must inspect tags() and dims() manually.\n");

  return m_data->dims[std::distance(m_data->tags.begin(),it)];
}

inline int FieldLayout::dim (const int idim) const {
  ekat::error::runtime_check(idim>=0 && idim<rank(), "Error! Index out of bounds.", -1);
  return m_data->dims[idim];
}

inline int FieldLayout::size () const {
  ekat::error::runtime_check(are_dimensions_set(), "Error! Field dimensions not yet set.\n",-1);
  int prod = rank()>0 ? 1 : 0;
  for (int idim=0; idim<rank(); ++idim) {
    prod *= m_data->dims[idim];
  }
  return prod;
}

inline FieldTag FieldLayout::tag (const int idim) const { 
  ekat::error::runtime_check(idim>=0 && idim<rank(), "Error! Index out of bounds.", -1);
  return m_data->tags[idim];
} 

inline bool FieldLayout::are_dimensions_set () const {
  return m_data->dims.size()>0;
}

inline bool operator== (const FieldLayout& fl1, const FieldLayout& fl2) {
  return fl1.rank()==fl2.rank() &&
         fl1.tags()==fl2.tags() &&
         fl1.dims()==fl2.dims();
}

} // namespace scream

#endif // SCREAM_FIELD_LAYOUT_HPP

