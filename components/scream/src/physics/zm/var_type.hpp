#ifndef VAR_TYPE_HPP
#define VAR_TYPE_HPP

#include <string>
#include <typeinfo>

std::string demangle(const char* name);

template <class T>
std::string var_type(const T& t) {

    return demangle(typeid(t).name());
}

#endif
