#ifndef SCREAM_FIELD_REPOSITORY_HPP
#define SCREAM_FIELD_REPOSITORY_HPP

#include "ekat/std_meta/ekat_std_utils.hpp"
#include "share/grid/grids_manager.hpp"
#include "share/field/field.hpp"
#include "share/field/field_request.hpp"
#include "share/util/map_key_iterator.hpp"
#include "share/scream_types.hpp"

#include "ekat/ekat_assert.hpp"
#include "ekat/util/ekat_string_utils.hpp"

#include <map>
#include <memory>
#include <set>

namespace scream
{

 /*
  *  A database for all the persistent fields needed in an atm time step
  *  We template a field repository over the field's (real) value type.
  *  This is enough to fully deduce the type of the stored views. All views
  *  are stored on the default device.
  *
  *  The fields are internally organized by name. Within each name,
  *  there can be multiple fields, which differ by their layout.
  *  For instance, we could have two version of 'temperature',
  *  one on the "physics" grid (tags: Column, Level), and one on
  *  the "dynamics" grid (tags: Element, GaussPoint, GaussPoint, Level).
  *
  *  When you query the repo for its 'size', you get the number
  *  of different *names*. So if the repo is storing the two
  *  versions of 'temperature' above and nothing else, its size
  *  will be 1. To get the number of different Field objects
  *  stored, you need to use the 'internal_size' method.
  *
  */

template<typename RealType>
class FieldRepository {
public:

  // Public types
  using RT               = typename std::remove_const<RealType>::type;
  using const_RT         = typename std::add_const<RT>::type;
  using field_type       = Field<RT>;
  using const_field_type = typename Field<RT>::const_field_type;
  using header_type      = typename field_type::header_type;
  using identifier_type  = typename field_type::identifier_type;
  using ci_string        = typename identifier_type::ci_string;
  using alias_map_type   = std::map<identifier_type,std::shared_ptr<field_type>>;
  using repo_type        = std::map<ci_string,alias_map_type>;
  using group_info_type  = FieldGroupInfo;
  using group_info_map   = std::map<ci_string,std::shared_ptr<group_info_type>>;

  // Constructor(s)
  FieldRepository ();

  // No copies, cause the internal database is not a shared_ptr.
  // NOTE: you can change this if you find that copies are needed/useful.
  FieldRepository (const FieldRepository&) = delete;
  FieldRepository& operator= (const FieldRepository&) = delete;

  // Change the state of the database
  void begin_registration ();
  void end_registration (const std::shared_ptr<const GridsManager>& gm = nullptr);
  void clean_up ();

  // Register a field in the repo.
  void register_field (const FieldRequest& req);

  // Set a GroupRequest in the repo. At registration_ends time, all requests will be
  // processed to try to meet demands. Notice that it might not be possible. E.g.,
  // the groups (A,B), (A,C), and (B,C) cannot be all available as a bundled field 
  void add_group_request (const GroupRequest& req);

  // Get information about the state of the repo
  int size () const { return m_fields.size(); }

  // Query for a particular field or group of fields
  bool has_field (const std::string& name) const;
  bool has_field (const identifier_type& identifier) const;
  field_type get_field (const identifier_type& identifier) const;
  field_type get_field(const std::string& name,const std::string& grid) const;

  // Get iterators to the keys (i.e., identifier_type) of all fields with a given name
  map_key_const_iterator<alias_map_type> cbegin (const std::string& name) const;
  map_key_const_iterator<alias_map_type> cend   (const std::string& name) const;

  // Get the names of the groups, with the names of all fields belonging to each group
  const group_info_map& get_groups_info () const { return m_field_groups; }

  // Note: when you request a group of fields, you must also specify the grid on which you need them.
  //       If you need the group G on grids foo and bar, then you must issue two separate requests to
  //       the field repo, for the group/grid pairs (G,foo) and (G,bar).
  // Note: it is *ASSUMED* that for each field in the group, there is only ONE field with such name
  //       on the requested grid.
  FieldGroup<RT> get_field_group (const std::string& name, const std::string& grid_name) const;
  FieldGroup<const_RT> get_const_field_group (const std::string& group_name, const std::string& grid_name) const;

  // Set the time stamp of all fields
  void init_fields_time_stamp (const util::TimeStamp& t0);

protected:

  std::shared_ptr<field_type> get_field_ptr(const std::string& name,const std::string& grid) const;
  std::shared_ptr<field_type> get_field_ptr(const identifier_type& id) const;

  // The state of the repository.
  RepoState   m_repo_state;

  // The actual repo.
  repo_type           m_fields;

  // The map group_name -> FieldGroupInfo
  group_info_map            m_field_groups;
  std::list<GroupRequest>   m_group_requests;

  // For each field, store the names of the groups it belongs to
  std::map<ci_string,std::set<ci_string>> m_field_to_groups;

  // These groups are allocated as a big bundled field, and each
  // field belonging to this group is created as a subview of that field
  // The map works m_bundled_groups[group_name] = bundled_field_name
  std::map<ci_string,ci_string>  m_bundled_groups;
};

// ============================== IMPLEMENTATION ============================= //

template<typename RealType>
FieldRepository<RealType>::
FieldRepository ()
 : m_repo_state (RepoState::Clean)
{
  m_bundled_groups["TRACERS"] = "Q";
  m_bundled_groups["TRACERS TENDENCY"] = "FQ";
  m_field_groups["TRACERS"] = std::make_shared<group_info_type>("TRACERS");
  m_field_groups["TRACERS TENDENCY"] = std::make_shared<group_info_type>("TRACERS TENDENCY");
}

template<typename RealType>
void FieldRepository<RealType>::
register_field (const FieldRequest& req) {
  // Sanity check
  EKAT_REQUIRE_MSG (m_repo_state==RepoState::Open,
      "Error! Repo state is not 'Open'. You should register fields between calls\n"
      "       to 'begin_registration()' and 'end_registration()'.\n");

  // Get the map of all fields with this name
  const auto& id = req.fid;
  auto& map = m_fields[id.name()];

  if (map.size()>0) {
    using ekat::units::to_string;
    // Make sure the new field id stores the same units as all the other ones.
    // TODO: this is the easiest way to ensure everyone uses the same units.
    //       However, in the future, we *may* allow different units, providing
    //       the users with conversion routines perhaps.
    const auto& id0 = map.begin()->first;
    EKAT_REQUIRE_MSG(id.get_units()==id0.get_units(),
        "Error! Request for field '" + id.name() + "' incompatible with already stored fields:\n"
        "        - input f units: " + to_string(id.get_units()) + "'\n"
        "        - stored f units: " + to_string(id0.get_units()) + "'\n"
        "      Please, check and make sure all atmosphere processes use the same units.\n");

    EKAT_REQUIRE_MSG (get_layout_type(id.get_layout().tags())==get_layout_type(id0.get_layout().tags()),
        "Error! Request for field '" + id.name() + "' incompatible with already stored fields:\n"
        "        - input f layout type:  " + to_string(get_layout_type(id.get_layout().tags())) + "\n"
        "        - stored f layout type: " + to_string(get_layout_type(id0.get_layout().tags())) + "\n"
        "      Please, check and make sure all atmosphere processes use compatible layouts.\n");

    // Silence compiler warnings in RELEASE mode
    (void) id0;
  }

  // Create (if needed), and get the field
  if (map.find(id)==map.end()) {
    map[id] = std::make_shared<field_type>(id);
  }
  auto& f = *map[id];

  // Make sure the field can accommodate the requested value type
  constexpr int real_size = sizeof(RealType);
  f.get_header().get_alloc_properties().request_allocation(real_size,req.pack_size);

  // Finally, add the field to the given groups
  // Note: we do *not* set the group info struct in the field header yet.
  //       we will do that when we end the registration phase.
  //       The reason is that when registration ends we will know *all* the groups
  //       that each field belongs to.
  for (const auto& group_name : req.groups) {
    // Get group (and init ptr, if necessary)
    auto& group = m_field_groups[group_name];
    if (group==nullptr) {
      group = std::make_shared<group_info_type>(group_name);
    }
    
    // Add the field name to the list of fields belonging to this group
    if (not ekat::contains(group->m_fields_names,id.name())) {
      group->m_fields_names.push_back(id.name());
    }
  }
}

template<typename RealType>
void FieldRepository<RealType>::
add_group_request (const GroupRequest& req) {
  EKAT_REQUIRE_MSG (m_repo_state==RepoState::Open,
      "Error! Repo state is not 'Open'. You should add group requests between calls\n"
      "       to 'begin_registration()' and 'end_registration()'.\n");

  // Add this req to the list of request.
  m_group_requests.push_back(req);
}

template<typename RealType>
bool FieldRepository<RealType>::
has_field (const std::string& name) const {
  auto it = m_fields.find(name);
  return it!=m_fields.end();
}

template<typename RealType>
bool FieldRepository<RealType>::
has_field (const identifier_type& identifier) const {
  auto it = m_fields.find(identifier.name());
  return it!=m_fields.end() && it->second.find(identifier)!=it->second.end();
}

template<typename RealType>
typename FieldRepository<RealType>::field_type
FieldRepository<RealType>::get_field (const identifier_type& id) const {
  EKAT_REQUIRE_MSG(m_repo_state==RepoState::Closed,
      "Error! Cannot get fields from the repo while registration has not yet completed.\n");
  auto ptr = get_field_ptr(id);
  EKAT_REQUIRE_MSG(ptr!=nullptr, "Error! Field not found.\n");
  return *ptr;
}

template<typename RealType>
typename FieldRepository<RealType>::field_type
FieldRepository<RealType>::get_field (const std::string& name, const std::string& grid) const {
  EKAT_REQUIRE_MSG(m_repo_state==RepoState::Closed,
      "Error! Cannot get fields from the repo while registration has not yet completed.\n");
  auto ptr = get_field_ptr(name,grid);
  EKAT_REQUIRE_MSG(ptr!=nullptr, "Error! Field " + name + " not found.\n");
  return *ptr;
}

template<typename RealType>
FieldGroup<typename FieldRepository<RealType>::RT>
FieldRepository<RealType>::
get_field_group (const std::string& group_name, const std::string& grid_name) const {
  EKAT_REQUIRE_MSG(m_repo_state==RepoState::Closed,
      "Error! Cannot get field groups from the repo while registration has not yet completed.\n");

  // Create an empty group
  FieldGroup<RT> group(group_name,grid_name);

  // Allow returning an empty group
  auto it = m_field_groups.find(group_name);
  if (it!=m_field_groups.end()) {
    group.m_info = it->second;;
    for (const auto& fname : group.m_info->m_fields_names) {
      auto f = get_field_ptr(fname,grid_name);
      group.m_fields[fname] = f;
      // Fetch the bundle fields (if bundled) just once
      if (group.m_info->m_bundled && group.m_bundle==nullptr) {
        auto p = f->get_header().get_parent().lock();
        EKAT_REQUIRE_MSG(p!=nullptr,
            "Error! A field belonging to a bundled field group is missing its 'parent'.\n");

        const auto& p_id = p->get_identifier();
        group.m_bundle = get_field_ptr(p_id);
      }
    }
  }

  return group;
}

template<typename RealType>
FieldGroup<typename FieldRepository<RealType>::const_RT>
FieldRepository<RealType>::
get_const_field_group (const std::string& group_name, const std::string& grid_name) const {
  EKAT_REQUIRE_MSG(m_repo_state==RepoState::Closed,
      "Error! Cannot get field groups from the repo while registration has not yet completed.\n");

  // Create an empty group
  FieldGroup<const_RT> group(group_name,grid_name);

  // Allow returning an empty group
  auto it = m_field_groups.find(group_name);
  if (it!=m_field_groups.end()) {
    group.m_info = it->second;;
    for (const auto& fname : group.m_info->m_fields_names) {
      auto f = get_field_ptr(fname,grid_name);
      auto cf = std::make_shared<const_field_type>(f->get_const());
      group.m_fields[fname] = cf;
      // Fetch the bundle fields (if bundled) just once
      if (group.m_info->m_bundled && group.m_bundle==nullptr) {
        auto p = f->get_header().get_parent().lock();
        EKAT_REQUIRE_MSG(p!=nullptr, "Error! Something is amiss with a bundled field group.\n");

        const auto& p_id = p->get_identifier();
        group.m_bundle = std::make_shared<const_field_type>(get_field_ptr(p_id)->get_const());
      }
    }
  }

  return group;
}

template<typename RealType>
void FieldRepository<RealType>::
init_fields_time_stamp (const util::TimeStamp& t0)
{
  EKAT_REQUIRE_MSG(m_repo_state==RepoState::Closed,
      "Error! Cannot set initial time stamp until registration has completed.\n");

  for (auto it : m_fields) {
    for (auto f : it.second) {
      f.second->get_header().get_tracking().update_time_stamp(t0);
    }
  }
}

template<typename RealType>
map_key_const_iterator<typename FieldRepository<RealType>::alias_map_type>
FieldRepository<RealType>::cbegin (const std::string& name) const
{
  EKAT_REQUIRE_MSG (m_fields.find(name)!=m_fields.end(),
      "Error! No field called '" << name << "' registered in the repo.\n");
  const auto& aliases = m_fields.at(name);
  return map_key_const_iterator<alias_map_type>(aliases.cbegin());
}

template<typename RealType>
map_key_const_iterator<typename FieldRepository<RealType>::alias_map_type>
FieldRepository<RealType>::cend (const std::string& name) const
{
  EKAT_REQUIRE_MSG (m_fields.find(name)!=m_fields.end(),
      "Error! No field called '" << name << "' registered in the repo.\n");
  const auto& aliases = m_fields.at(name);
  return map_key_const_iterator<alias_map_type>(aliases.cend());
}

template<typename RealType>
void FieldRepository<RealType>::begin_registration ()
{
  // Update the state of the repo
  m_repo_state = RepoState::Open;
}

template<typename RealType>
void FieldRepository<RealType>::
end_registration (const std::shared_ptr<const GridsManager>& gm)
{
  // Ugly complexity, but we won't have that many fields, and this is all setup costs.
  auto contains = [] (const std::list<ci_string>& super,
                      const std::list<ci_string>& sub) -> bool {
    for (const auto& s : sub) {
      if (not ekat::contains(super,s)) {
        return false;
      }
    }
    // All items in sub are also in super.
    return true;
  };
  auto intersects = [] (const std::list<ci_string>& l1,
                        const std::list<ci_string>& l2) -> bool {
    for (const auto& s : l1) {
      if (ekat::contains(l2,s)) {
        return true;
      }
    }
    // All items in l1 are not in l2.
    return false;
  };

  // Before allocating fields, take a moment to gather groups information.
  // In particular, if a group is requested bundled, then we need to fist
  // allocate the group, and then subview the individual fields.
  // Things get complicated if there are two groups that overlap. In particular,
  // the hard case is if group G1 and G2 share some fields, but not all.
  // In this case, we need to check if it's possible to bundle both groups or not.
  // If not possible, and both were required 'bundled', we crap out.
  // With 2 groups, it should always be possible, but with 3, it might not be.
  // E.g., G1=(A,B), G2=(B,C), G3=(A,C) cannot all be bundled at the same time,
  // but any two of them can.

  // Checks on group requests
  for (const auto& req : m_group_requests) {
    // Check this group exists
    EKAT_REQUIRE_MSG (m_field_groups.find(req.name)!=m_field_groups.end(),
        "Error! Request for non-existent group '" + req.name + "'.\n");

    if (req.superset_group!="") {
      // Check the superset group exists
      auto super = m_field_groups.find(req.superset_group);
      EKAT_REQUIRE_MSG (not super==m_field_groups.end(),
          "Error! Request for group '" + req.name + "' is invalid.\n"
          "       Non-existing superset group '" + req.superset_group + "'.\n");

      // Check the fields to exclude from superset are indeed in the superset
      for (const auto& ex : req.exclude_superset_fields) {
        EKAT_REQUIRE_MSG (ekat::contains(super->second.m_fields_names,ex),
            "Error! Request for group '" + req.name + "' is invalid.\n"
            "       Field '" + ex + "' not found in superset group '" + req.superset_group + "'.\n");
      }
    }
  }

  // For each group, pick the strongest request in terms of bundling.
  std::map<ci_string,GroupBundling> bundling;
  std::map<ci_string,std::set<ci_string>> grids;
  for (const auto& req : m_group_requests) {
    bundling[req.name] = std::max(bundling[req.name],req.bundling);
    grids[req.name].insert(req.grid);
  }

  // If a group is to be bundled, all fields must be "compatible".
  // E.g., we can't bundle a 2d field and a 3d one.
  auto get_field_dim = [] (const FieldLayout& layout) -> int {
    auto type = get_layout_type(layout.tags());
    if (type==LayoutType::Scalar2D || type==LayoutType::Vector2D) {
      return 2;
    } else if (type==LayoutType::Scalar3D || type==LayoutType::Vector3D) {
      return 3;
    } else {
      // Only support 2d and 3d, scalar or vector
      EKAT_ERROR_MSG ("Error! Unsupported layout '" + e2str(type) +"'.\n");
    }
  };
  for (auto it : bundling) {
    if (it.second!=NotNeeded) {
      int dim = -1;
      const auto& names = m_field_groups.at(it.first);
      for (const auto& n : names) {
        // Pick any field with this name (doesn't matter the grid);
        const auto& f = *m_fields.at(n).begin();
        const auto fdim = get_field_dim(f.get_header().get_identifier().get_layout());
        if (dim==-1) {
          dim = fdim;
        } else {
          EKAT_REQUIRE_MSG (fdim==dim,
              "Error! Request for group '" + it.first + "' to be bundled cannot be fulfilled.\n"
              "       Group contains both 3d and 2d fields.\n");
        }
      }
    }
  }

  // This is the *hard* part. We need to figure out if 2+ groups overlap.
  // If they do *and* they all need to be bundled, then we need to check
  // if we can rearrange the field_names in all the groups, in a way that
  // accommodates all the groups.
  // We support only three cases:
  //  - G1 subset of G2. Easy: put all G1 first (or last) in G2
  //  - G1 intersect G2, but no subset relationship. Here, we create
  //    a "support" group G3=G1+G2. Then, in G3 we put G1-G2 first,
  //    then the intersection, and finally G2-G1.
  //  - G1 contains G2...Gn, but the Gi's intersect only pairwise,
  //    that is, G1 intersects G2, G3 intersects G4,...
  //    Here, do the step above processing the groups 2 at a time.
  // There are other scenarios that can be accommodate, but we defer their
  // implementation to when it is actually needed.
  // CAVEAT: the 'tracers' group adds a bit of complexity, since Homme
  // (and potentially other parts of EAM) assume that water-vapor is
  // the 1st tracer in the bundled array Q. So we need to hard-code
  // Q=(qv,...).

  // Establish parent/sibling relationships. G1 is a parent of G2 if it
  // contains *at least* all fields in G2. G1 is a sibling of G2 if they
  // have common fields, but none is parent of the other.
  // TODO: there might be legit cases where G1==G2 (they are aliases), but
  //       for now we just don't allow it.
  // Note: if a group has bundling=NotNeeded, we don't bother finding its
  //       parent(s), nor its siblings. And we don't count it as a sibling
  //       for other groups.
  std::map<ci_string,std::list<ci_string>> parents, siblings;
  for (const auto& it1 : m_field_groups) {
    for (const auto& it2 : m_field_groups) {
      if (it1.first==it2.first) {
        // skip self-comparison
        continue;
      }
      if (bundling[it1.first]==NotNeeded) {
        // This group is fine with whatever order the fields are
      }

      const auto& g1 = it1.second->m_fields_names;
      const auto& g2 = it2.second->m_fields_names;
      if (contains(g2,g1)) {
        EKAT_REQUIRE_MSG(not contains(g2,g1),
            "Error! We do not allow a group to alias another group.\n");
        parents[it1.first].push_back(it2.first);
      } else if (bundling[it2.first]!=NotNeeded and
                 intersects(g1,g2) and not contains(g2,g1)) {
        // Store sibling only if both might need bundling
        siblings[it1.first].push_back(it2.first);
      }
    }
  }

  // For ease of implementation, we require at most one parent and one
  // sibling per group. That means that for each G1 there's at most
  // one G2 that contains G1, and at most one G3 that intersects G1.
  for (const auto& it : parents) {
    if (it.second.size()>=0) {
      EKAT_REQUIRE_MSG (it.second.size()<=1,
          "Error! We do not yet support nested field groups.\n");

      // If this group needs bundling, its parent better be bundled too
      bundling[it.second.front()] = std::max(bundling[it.second.front()],
                                             bundling[it.first]);
    }

    // Make sure we intersect with at most one other group
    // If there are N>1 siblings, discard N-1 siblings for which
    // bundling is optional, if possible. If we still end up with
    // 2+ siblings, crap out.
    auto& sib = siblings[it.first];
    if (sib.size()>1) {
      auto bundling_optional = [&](const ci_string& name) {
        return bundling[name]!=Needed;
      };
      auto pos = std::find(sib.begin(),sib.end(),bundling_optional);
      while (pos!=sib.end() && sib.size()>1) {
        sib.erase(pos);
        pos = std::find(sib.begin(),sib.end(),bundling_optional);
      }
      EKAT_REQUIRE_MSG (it.second.size()<=1,
          "Error! We do not yet support intersection of 3+ bundled groups of fields.\n");
    }
  }

  // If two groups intersect, and do not have a parent, we create
  // a new "support" group, G1+G2. Notice that, since we only allow *one*
  // intersection, G1+G2 cannot intersect any other group. Also, since
  // we only allow one parent and one sibling, either G1 and G2 have the same parent,
  // or neither of them has a parent.
  for (const auto& it : siblings) {
    const auto& g1_name = it.first;

    if (it.second.size()>0) {
      const auto& g2_name = it.second.front();
      auto& p1 = parents[g1_name];
      auto& p2 = parents[g2_name];

      // This check should be implied by the previous ones, but better be safe
      EKAT_REQUIRE_MSG (p1.size()==p2.size(),
          "Error! I found two groups that intersect, but one has a parent, and the other doesn't.\n");

      if (p1.size()==0) {
        // Create the FieldGroupInfo
        ci_string g1_g2_name = "union_" + g1_name + "_" + g2_name;
        auto& g1_g2 = m_field_groups[g1_g2_name] = std::make_shared<FieldGroupInfo>(g1_g2_name);
        const auto& g1 = m_field_groups[g1_name];
        const auto& g2 = m_field_groups[g2_name];
        g1_g2->m_fields_names = g1->m_fields_names;
        g1_g2->m_fields_names.insert(g1_g2->m_fields_names.end(),g2->m_fields_names.begin(),g2->m_fields_names.end());
        g1_g2->m_fields_names.sort();
        g1_g2->m_fields_names.unique();

        // Add g1_g2 as parent of g1 and g2
        p1.push_back(g1_g2_name);
        p2.push_back(g1_g2_name);

        // Set specs for g1_g2
        bundling[g1_g2_name] = std::max(bundling[g1_name],bundling[g2_name]);
        grids[g1_g2_name] = grids[g1_name];
        grids[g1_g2_name].insert(grids[g2_name].begin(),grids[g2_name].end());
      }
    }
  }

  // Go over the groups, and reorder fields
  std::map<ci_string,std::list<ci_string>> sorted_names;
  std::map<ci_string,bool> done;
  for (auto& it : m_field_groups) {
    const auto& name = it.first;

    // If no bundling needed, we don't care how fields are ordered
    if (bundling[name]==NotNeeded) {
      continue;
    }

    // We do not intersect with another group, so we don't care how fields are ordered
    // Also, skip this field if we already processed it when we processed its sibling
    if (siblings[name].size()==0 || done[name]) {
      continue;
    }

    // This group (G1) is overlapping with G2. Reorder our fields
    // to get G1=(f1,f12), G2=(f12,f2), and add (f1,f12,f2) to
    // our parent list of fields.
    auto& s = m_field_groups[siblings[name]];
    std::list<ci_string>& g1 = it.second->m_fields_names;
    std::list<ci_string>& g2 = s->m_fields_names;

    auto in_g1 = [&](const ci_string& s) {
      return ekat::contains(g1,s);
    };
    auto not_in_g2 = [&](const ci_string& s) {
      return ekat::contains(g2,s);
    };
    std::partition(g1.begin(),g1.end(),not_in_g2);
    std::partition(g2.begin(),g2.end(),in_g1);
    sorted_names[name] = g1;
    sorted_names[siblings[name]] = g2;

    // Now, add g1+g2 (removing duplicates) in the list of fields of our parent
    const auto& p = parents[name].front();
    sorted_names[p] = g1;
    for (const auto& n : g2) {
      if (not ekat::contains(sorted_names[p],n)) {
        sorted_names[p].push_back(n);
      }
    }

    done[it.first] = true;
    done[siblings[name]] = true;
  }

  // Parent groups might still have some fields not inserted in the fields_names lists.
  // That is the case if, say, G1=(C,B), G2=(A,D,C), G3=(A,B,C,D,E). After the above
  // loop, sorted_names[G3]=(C,B,A,D), so we need to add all remaining fields at the end.
  for (auto& it : m_field_groups) {
    auto& fnames = it.second->m_fields_names;
    auto& sorted = sorted_names[it.first];
    if (!done[it.first]) {
      for (const auto& n : fnames) {
        if (!ekat::contains(sorted,n)) {
          sorted.push_back(n);
        }
      }
    }

    // Reset the m_fields_names to store the sorted names
    fnames = sorted;
  }

  // Count fields in each group, and list the grids they are needed on
  std::map<ci_string,int> sizes;
  for (const auto& it : m_field_groups) {
    const auto& names = it.second->m_fields_names;
    int& s = sizes[it.first];
    s = 0;
    for (const auto& n : names) {
      const auto& f = *m_fields.at(n).begin();

      const auto& layout = f.get_header().get_identifier().get_layout();
      const auto lt = get_layout_type(layout.tags());
      if (lt==LayoutType::Vector2D || lt==LayoutType::Vector3D) {
        s += layout.dim(FieldTag::Component);
      } else {
        s += 1;
      }
    }
  }

  // FINALLY we have all the names for all the groups, in an order that satisfy bundle requests.
  // First, allocate groups that have no parent
  for (const auto& it : m_field_groups) {
    const auto& name = it.first;
    if (parents[it.first].size()>0) {
      continue;
    }
    ///////////////////////////
    for (const auto& grid : grids[name]) {
      if (bundling[name]) {
        // Create field id for Q
        auto layout = gm->get_grid(gn)->get_3d_vector_layout(true,VAR,nt);
        FieldIdentifier fid_Q (Q_name,layout,q_units,gn);

        // Create Q field
        register_field(fid_Q);
        auto Q = get_field_ptr(fid_Q);

        // Scan all tracers on this grid, get their alloc prop,
        // and make sure Q can accommodate all of them
        auto& Q_ap = Q->get_header().get_alloc_properties();
        for (const auto& fn : tr_gr.m_fields_names) {
          auto q = get_field_ptr (fn, gn);
          if (q!=nullptr) {
            Q_ap.request_allocation(q->get_header().get_alloc_properties());
          }
        }

        // Allocate
      } else {
      }
    }
  }

  // If there are tracers, we need gm to be valid
  EKAT_REQUIRE_MSG (nt==0 || gm!=nullptr,
      "Error! Tracers allocation requires a valid grid manager pointer.\n");

  const std::string& Q_name = m_bundled_groups.at("TRACERS");
  const std::string& FQ_name = m_bundled_groups.at("TRACERS TENDENCY");
  constexpr auto VAR = ShortFieldTagsNames::VAR;
  const auto q_units = ekat::units::Units::nondimensional();

  // Create "bundled" tracers and tracers forcing
  for (const auto& gn : tr_grids) {
    // Create field id for Q
    auto layout = gm->get_grid(gn)->get_3d_vector_layout(true,VAR,nt);
    FieldIdentifier fid_Q (Q_name,layout,q_units,gn);

    // Create Q field
    register_field(fid_Q);
    auto Q = get_field_ptr(fid_Q);

    // Scan all tracers on this grid, get their alloc prop,
    // and make sure Q can accommodate all of them
    auto& Q_ap = Q->get_header().get_alloc_properties();
    for (const auto& fn : tr_gr.m_fields_names) {
      auto q = get_field_ptr (fn, gn);
      if (q!=nullptr) {
        Q_ap.request_allocation(q->get_header().get_alloc_properties());
      }
    }

    // Allocate
    Q->allocate_view();
  }

  // Helper lambdas to detect if a field name corresponds to a tracer or
  // tracer tendency.
  auto is_tracer = [&](const std::string& name) -> bool {
    return ekat::contains(tr_gr.m_fields_names,name);
  };

  // Proceed to allocate other fields, and subview tracers
  for (auto& it : m_fields) {
    const auto& fname = it.first;
    if (fname==Q_name || fname==FQ_name) {
      // We already allocated Q and FQ, so skip it.
      continue;
    }

    for (auto& it_f : it.second) {
      auto f = it_f.second;
      if (is_tracer(fname)) {
        // A tracer must be a subview of the big Q field.
        auto pos = ekat::find(tr_gr.m_fields_names,fname);
        const int iq = std::distance(tr_gr.m_fields_names.begin(),pos);
        EKAT_REQUIRE_MSG (iq>=0 && iq<tr_gr.size(),
            "Error! Field '" << fname << "' is a tracer, but could not locate it in the tracer group.\n");

        const auto& fh = f->get_header();
        const auto  Q = get_field_ptr(Q_name,fh.get_identifier().get_grid_name());
        const auto& Q_tags = Q->get_header().get_identifier().get_layout().tags();
        // Note: as of 02/2021, idim should *always* be 1, but we store it just in case,
        //       to avoid bugs in the future.
        const int idim = std::distance(Q_tags.begin(),ekat::find(Q_tags,VAR));
        auto q = Q->subfield(fname,fh.get_identifier().get_units(),idim,iq);

        // Either this is the first tracer we set in the group (m_subview_dim still -1),
        // or idim should match what was already in the group info.
        EKAT_REQUIRE_MSG (tr_gr.m_subview_dim==-1 || tr_gr.m_subview_dim==idim,
            "Error! Something is amiss with the creation of tracers subviews.\n");
        tr_gr.m_subview_idx[fname] = iq;
        tr_gr.m_subview_dim = idim;
        tr_gr.m_bundled = true;

        // Overwrite f with q.
        *f = q;
      } else {
        // A completely independent field. Allocate it.
        f->allocate_view();
      }
    }
  }

  // Update the tracking of all fields with the group info of all
  // the groups they belong to
  for (const auto& fgi_it : m_field_groups) {
    auto fgi = fgi_it.second;
    // Get fields in this group
    const auto& fnames = fgi->m_fields_names;
    for (const auto& fn : fnames) {
      // Get all fields with this name
      auto& map = m_fields.at(fn);
      for (auto& it : map) {
        // Update the field tracking
        it.second->get_header().get_tracking().add_to_group(fgi);
      }
    }
  }

  // Prohibit further registration of fields
  m_repo_state = RepoState::Closed;
}

template<typename RealType>
void FieldRepository<RealType>::clean_up() {
  m_fields.clear();
  m_repo_state = RepoState::Clean;
}

template<typename RealType>
std::shared_ptr<typename FieldRepository<RealType>::field_type>
FieldRepository<RealType>::get_field_ptr (const identifier_type& id) const {
  if (!has_field(id)) {
    return nullptr;
  }
  return m_fields.at(id.name()).at(id);
}

template<typename RealType>
std::shared_ptr<typename FieldRepository<RealType>::field_type>
FieldRepository<RealType>::get_field_ptr (const std::string& name, const std::string& grid) const {

  if (!has_field(name)) {
    return nullptr;
  }

  // Keep track of the number of fields found for this name/grid combo
  std::vector<identifier_type> f_matches;

  //  Search subset of field repo for matching gridname.
  for (const auto& it : m_fields.at(name)) {
    if ( it.first.get_grid_name()==grid) {
      f_matches.push_back(it.first);
    }
  }

  // Check to make sure a) the field was found on this grid, and b) only once.
  EKAT_REQUIRE_MSG(f_matches.size()==1, "Error! get_field: " + name + " found " + std::to_string(f_matches.size()) + " matches on grid " + grid + ".\n");

  // Use this field id to grab field itself.
  return get_field_ptr(f_matches[0]);
}

} // namespace scream

#endif // SCREAM_FIELD_REPOSITORY_HPP
