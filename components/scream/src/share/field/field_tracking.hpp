#ifndef SCREAM_FIELD_TRACKING_HPP
#define SCREAM_FIELD_TRACKING_HPP

#include "share/field/field_group.hpp"
#include "share/scream_types.hpp"
#include "share/util/scream_time_stamp.hpp"
#include "share/util/scream_family_tracking.hpp"

#include "ekat/util/ekat_string_utils.hpp"
#include "ekat/std_meta/ekat_std_utils.hpp"
#include "ekat/ekat_assert.hpp"

#include <memory>   // For std::weak_ptr
#include <string>
#include <list>
#include <set>

namespace scream {

// Forward declarations
class AtmosphereProcess;
class FieldHeader;

class FieldTracking : public FamilyTracking<FieldTracking> {
public:

  using TimeStamp         = util::TimeStamp;
  using ci_string         = ekat::CaseInsensitiveString;
  using atm_proc_ptr_type = std::weak_ptr<AtmosphereProcess>;
  using atm_proc_set_type = ekat::WeakPtrSet<AtmosphereProcess>;

  FieldTracking() = delete;
  explicit FieldTracking(const std::string& name);
  FieldTracking(const FieldTracking&) = default;

  // No assignment, to prevent tampering with tracking (e.g., rewinding time stamps)
  FieldTracking& operator=(const FieldTracking&) = delete;

  // ----- Getters ----- //

  // The time stamp of the field. This can be used to check when it was last updated.
  // Please, notice this is not the OS time stamp (see time_stamp.hpp for details).
  const TimeStamp& get_time_stamp () const { return m_time_stamp; }

  // List of providers/customers for this field
  const atm_proc_set_type& get_providers () const { return m_providers; }
  const atm_proc_set_type& get_customers () const { return m_customers; }

  // List of field groups that this field belongs to
  const ekat::WeakPtrSet<const FieldGroupInfo>& get_groups_info () const { return m_groups; }

  // ----- Setters ----- //

  // Add to the list of providers/customers
  void add_provider (const std::weak_ptr<AtmosphereProcess>& provider);
  void add_customer (const std::weak_ptr<AtmosphereProcess>& customer);

  // Add the field to a given group
  void add_to_group (const std::shared_ptr<const FieldGroupInfo>& group);

  // Set the time stamp for this field. This can only be called once, due to TimeStamp implementation.
  // NOTE: if the field has 'children' (see FamilyTracking), their ts will be updated too.
  //       However, if the field has a 'parent' (see FamilyTracking), the parent's ts will not be updated.
  void update_time_stamp (const TimeStamp& ts);

  const std::string& name () const { return m_name; }

protected:

  // We keep the field name just to make debugging messages more helpful
  std::string m_name;

  // Tracking the updates of the field
  TimeStamp         m_time_stamp;

  // These are to be used to track the order in which providers update the field at each time step.
  // One can use this information to track when a field gets updated during a timestep. It can be
  // particularly useful in the case of parallel schedules.
  std::set<std::string>   m_last_timestep_providers;
  std::set<std::string>   m_curr_timestep_providers;

  // List of provider/customer processes. A provider is an atm process that computes/updates the field.
  // A customer is an atm process that uses the field just as an input.
  // NOTE: do NOT use shared_ptr, since you would create circular references.
  atm_proc_set_type       m_providers;
  atm_proc_set_type       m_customers;

  // Groups are used to bundle together fields, so that a process can request all of them
  // without knowing/listing all their names. For instance, the dynamics process needs to
  // get all tracers, which need to be advected. However, dyamics has no idea of what are
  // the tracers names, and neither should it care. Groups can come to rescue here, allowing
  // dynamics to request all fields that have been marked as 'tracers'.
  ekat::WeakPtrSet<const FieldGroupInfo>    m_groups;
};

// Use this free function to exploit features of enable_from_this
template<typename... Args>
inline std::shared_ptr<FieldTracking>
create_tracking(const Args&... args) {
  auto ptr = std::make_shared<FieldTracking>(args...);
  ptr->setSelfPointer(ptr);
  return ptr;
}


} // namespace scream

#endif // SCREAM_FIELD_TRACKING_HPP
