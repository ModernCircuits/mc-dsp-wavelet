#pragma once

#include "lt/preprocessor.hpp"

#if __has_include(<bit>)
    #include <bit>
#endif

#if defined(__cpp_lib_bit_cast)
namespace lt
{
using std::bit_cast;
}
#else
namespace lt
{

}  // namespace lt
#endif
