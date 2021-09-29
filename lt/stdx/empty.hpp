#pragma once

#include <initializer_list>
#include <iterator>
#include <type_traits>

#include "lt/preprocessor.hpp"

#if defined(__cpp_lib_nonmember_container_access)
namespace lt {
using std::empty;
}
#else
namespace lt {

template <typename C>
LT_NODISCARD constexpr auto empty(const C& c) -> decltype(c.empty())
{
    return c.empty();
}

template <typename T, std::size_t N>
LT_NODISCARD constexpr auto empty(const T (&array)[N]) noexcept -> bool
{
    (void)array;
    return false;
}

template <typename E>
LT_NODISCARD constexpr auto empty(std::initializer_list<E> il) noexcept -> bool
{
    return il.size() == 0;
}

} // namespace lt
#endif