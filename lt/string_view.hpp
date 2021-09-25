#pragma once

#if defined(__cpp_lib_string_view)
#include <string_view>
namespace lt {
using std::string_view;
} // namespace lt
#else
#include "boost/utility/string_view.hpp"
namespace lt {
using boost::string_view;
} // namespace lt
#endif