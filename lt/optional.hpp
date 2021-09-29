#pragma once

#if __has_include(<optional>)
#include <optional>
#endif

#if defined(__cpp_lib_optional)
namespace lt {
using std::optional;
} // namespace lt
#else
#include "boost/optional.hpp"
namespace lt {
using boost::optional;
} // namespace lt
#endif