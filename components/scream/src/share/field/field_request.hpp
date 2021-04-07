#ifndef SCREAM_FIELD_REQUEST_HPP
#define SCREAM_FIELD_REQUEST_HPP

#include "share/field//field_identifier.hpp"

#include <ekat/util/ekat_string_utils.hpp>
#include <ekat/util/ekat_units.hpp>
#include <ekat/ekat_assert.hpp>

#include <list>

namespace scream {

/*
 * A struct used to request a field.
 *
 * FieldRequest is a lightweight struct that an Atmosphere Process (AP) can
 * expose to signal the driver that it needs a certaing field.
 * The request must include the FieldIdentifier of the field, the desired
 * pack size for this field, and a list of groups to which the field is
 * meant to be associated with.
 */

struct FieldRequest {
  using ci_string  = ekat::CaseInsensitiveString;

  // Main ctor
  FieldRequest (const FieldIdentifier& fid_, const int ps, const std::list<ci_string>& groups_ = {})
   : fid(fid_), pack_size(ps), groups(groups_)
  {
    EKAT_REQUIRE_MSG(pack_size>=1, "Error! Invalid pack size request.\n");
  }

  // Other convenience ctors (they all call the 1st)
  FieldRequest (const FieldIdentifier& fid_, const std::list<ci_string>& groups_ = {})
   : FieldRequest(fid_,1,groups_)
  { /* Nothing to do here */ }

  FieldRequest (const FieldIdentifier& fid_, const int ps, const std::string& group)
   : FieldRequest(fid_,ps,std::list<ci_string>{group})
  { /* Nothing to do here */ }

  FieldRequest (const std::string& name, const FieldLayout& layout, const ekat::units::Units& u, const std::string& grid, const int ps, const std::list<ci_string>& groups_)
   : FieldRequest(FieldIdentifier(name,layout,u,grid),ps,groups_)
  { /* Nothing to do here */ }

  FieldRequest (const std::string& name, const FieldLayout& layout, const ekat::units::Units& u, const std::string& grid, const int ps, const std::string& group)
   : FieldRequest(FieldIdentifier(name,layout,u,grid),ps,std::list<ci_string>{group})
  { /* Nothing to do here */ }

  FieldIdentifier fid;
  int pack_size;
  std::list<ci_string> groups;

  // TODO: if/when we do mixed precision, add a member
  //   int real_size;
  // (or maybe an instance of enum Precision { single, double };
};

/*
 * A struct used to request a group of fields.
 *
 * A GroupRequest is a lightweight struct that an Atmosphere Process (AP) can
 * expose to signal the driver that it needs a certaing group of fields,
 * without the need of knowing how many fields are in the group, or how they are called.
 * A typical example is an AP that needs to advect tracers (like Dynamics does):
 * it treats tracers agnostically, and does not really care how many there are.
 * So the AP exposes this need as a GroupRequest. Later, it will be provided
 * with a FieldGroup, which allows to access all the fields in the group
 * individually, and, if the allocation permits it, as a single N+1 dimensional
 * field. For more details about the FieldGroup struct, see field_group.hpp.
 */

enum GroupBundling : int {
  // Values arranged in increasing order of demand
  NotNeeded = 0,
  Preferred = 1,
  Needed    = 2
};

struct GroupRequest {
  using ci_string  = ekat::CaseInsensitiveString;

  // Main ctor
  GroupRequest (const std::string& name_, const std::string& grid_,
                const int ps, const GroupBundling bundling_,
                const std::string& superset_group_,
                const std::list<std::string>& exclude_superset_fields_)
   : name(name_), grid(grid_), pack_size(ps), bundling(bundling_)
  {
    superset_group = superset_group_;
    EKAT_REQUIRE_MSG(pack_size>=1, "Error! Invalid pack size request.\n");
    for (const auto& n : exclude_superset_fields_) {
      exclude_superset_fields.push_back(n);
    }
  }

  // Other convenience ctors (they all call the 1st)
  GroupRequest (const std::string& name_, const std::string& grid_, const int ps, const GroupBundling bundling_)
   : GroupRequest(name_,grid_,ps,bundling_,"",{})
  { /* Nothing to do here */ }

  GroupRequest (const std::string& name_, const std::string& grid_)
   : GroupRequest(name_,grid_,1,NotNeeded)
  { /* Nothing to do here */ }

  GroupRequest (const std::string& name_, const std::string& grid_,
                const std::string& superset_group_,
                const std::list<std::string>& exclude_superset_fields_)
   : GroupRequest(name_,grid_,1,NotNeeded,superset_group_,exclude_superset_fields_)
  { /* Nothing to do here */ }

  // Group name
  ci_string name;
  // Grid name
  ci_string grid;
  // Request an allocation that can accomodate a value type like Pack<Real,pack_size>
  int       pack_size;
  // Whether a bundled field version of this group is needed, preferred, or not needed.
  // If needed, and the Repo cannot accommodate the request, an error will be thrown
  GroupBundling  bundling;

  // Allow to specify a group as a subset of another, excluding a given list of fields
  ci_string superset_group;
  std::list<ci_string> exclude_superset_fields;
};

} // namespace scream

#endif // SCREAM_FIELD_REQUEST_HPP
