#pragma once

#if defined(__cpp_lib_optional)
#include <optional>
namespace lt {
using std::optional;
} // namespace lt
#else
#include "boost/optional.hpp"
namespace lt {
using boost::optional;
} // namespace lt
#endif