#pragma once

#if __has_include(<filesystem>)
#include <filesystem>
#endif

#if defined(__cpp_lib_filesystem)
namespace lt {
namespace filesystem = std::filesystem;
} // namespace lt
#else
#include "boost/filesystem.hpp"
namespace lt {
namespace filesystem = boost::filesystem;
} // namespace lt
#endif