
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
LT_NODISCARD constexpr auto size(C const& c) -> decltype(c.size())
{
    return c.size();
}
template <typename T, std::size_t N>
LT_NODISCARD constexpr auto size(T const (&array)[N]) noexcept -> std::size_t
{
    (void)array;
    return N;
}
} // namespace lt
#endif