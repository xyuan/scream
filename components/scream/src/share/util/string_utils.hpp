#ifndef SCREAM_STRING_UTILS_HPP
#define SCREAM_STRING_UTILS_HPP

#include <string>
#include <sstream>
#include <algorithm>

namespace scream {
namespace util {

// Small utility that cats a space and an integer to an input string.
inline std::string strint (const std::string& s, const int i) {
  std::stringstream ss;
  ss << s << " " << i;
  return ss.str();
}

// A no-overhead class that inherits from std::string, which we only
// use to trigger different behavior in the ==,!=,<,<= operators.
class CaseInsensitiveString final : public std::string
{
public:
  template<typename... Args>
  CaseInsensitiveString (Args... args)
   : std::string(args...)
  {}

  virtual ~CaseInsensitiveString () = default;
};

// Case-insensitive comparison functions
bool caseInsensitiveEqualString (const std::string& s1, const std::string& s2);
bool caseInsensitiveLessEqualString (const std::string& s1, const std::string& s2);
bool caseInsensitiveLessEqualString (const std::string& s1, const std::string& s2);

inline bool caseInsensitiveEqualString (const std::string& s1, const std::string& s2) {
  auto charComp = [](const char c1, const char c2)->bool{
    return c1==c2 || std::toupper(c1)==std::toupper(c2);
  };
  return s1.size()==s2.size() &&
         std::equal(s1.begin(),s1.end(),s2.begin(),charComp);
}

inline bool caseInsensitiveLessString (const std::string& s1, const std::string& s2) {
  if (s1.size()<=s2.size()) {
    auto charComp = [](const char c1, const char c2)->bool{
      return std::toupper(c1)<std::toupper(c2);
    };

    return std::equal(s1.begin(),s1.end(),s2.begin(),charComp);
  } else {
    return !caseInsensitiveLessEqualString(s2,s1);
  }
}

inline bool caseInsensitiveLessEqualString (const std::string& s1, const std::string& s2) {
  if (s1.size()<=s2.size()) {
    auto charComp = [](const char c1, const char c2)->bool{
      return std::toupper(c1)<=std::toupper(c2);
    };
    return std::equal(s1.begin(),s1.end(),s2.begin(),charComp);
  } else {
    return !caseInsensitiveLessString(s2,s1);
  }
}

// Overloads of comparison operators, which use the routines above if at least one
// of the two inputs is indeed a CaseInsensitiveString
template<typename S1, typename S2>
typename std::enable_if<
      std::is_same<S1,CaseInsensitiveString>::value ||
      std::is_same<S2,CaseInsensitiveString>::value,
      bool>::type
operator== (const S1& s1, const S2& s2) {
  return caseInsensitiveEqualString(s1,s2);
}

template<typename S1, typename S2>
typename std::enable_if<
      std::is_same<S1,CaseInsensitiveString>::value ||
      std::is_same<S2,CaseInsensitiveString>::value,
      bool>::type
operator!= (const S1& s1, const S2& s2) {
  return ! (s1==s2);
}

template<typename S1, typename S2>
typename std::enable_if<
      std::is_same<S1,CaseInsensitiveString>::value ||
      std::is_same<S2,CaseInsensitiveString>::value,
      bool>::type
operator< (const S1& s1, const S2& s2) {
  return caseInsensitiveLessString(s1,s2);
}

template<typename S1, typename S2>
typename std::enable_if<
      std::is_same<S1,CaseInsensitiveString>::value ||
      std::is_same<S2,CaseInsensitiveString>::value,
      bool>::type
operator<= (const S1& s1, const S2& s2) {
  return caseInsensitiveLessEqualString(s1,s2);
}

} // namespace util
} // namespace scream

#endif // SCREAM_STRING_UTILS_HPP