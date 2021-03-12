#include "field_layout.hpp"

namespace scream
{

FieldLayout::FieldLayout ()
{
  m_data = std::make_shared<Data>();
}

FieldLayout::FieldLayout (const std::vector<FieldTag>& tags)
 : FieldLayout ()
{
  m_data->tags = tags;
}

FieldLayout::FieldLayout (const std::vector<FieldTag>& tags,
                          const std::vector<int>& dims)
 : FieldLayout (tags)
{
  set_dimensions(dims);
}

void FieldLayout::set_dimensions (const std::vector<int>& dims) {
  // Check, then set dims
  EKAT_REQUIRE_MSG (!are_dimensions_set(),
      "Error! You cannot reset dimensions once they have been set.\n");
  EKAT_REQUIRE_MSG(dims.size()==m_data->tags.size(),
      "Error! Input dimensions vector not properly sized.");

  m_data->dims.resize(m_data->tags.size());
  for (int idim=0; idim<rank(); ++idim) {
    m_data->dims[idim] = dims[idim];
  }
}

} // namespace scream
