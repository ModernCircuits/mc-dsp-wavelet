
#pragma once

#include <iterator>
#include <type_traits>

#include "lt/preprocessor.hpp"

#if defined(__cpp_lib_nonmember_container_access)
namespace lt {
using std::size;
}
#else
namespace lt {
template <typename C>
LT_NODISCARD constexpr auto size(const C& c) -> decltype(c.size())
{
    return c.size();
}
template <typename T, std::size_t N>
LT_NODISCARD constexpr auto size(const T (&array)[N]) noexcept -> std::size_t
{
    return N;
}
} // namespace lt
#endif