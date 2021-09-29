#pragma once

#include <initializer_list>
#include <iterator>
#include <type_traits>

#include "lt/preprocessor.hpp"

#if defined(__cpp_lib_nonmember_container_access)
namespace lt {
using std::data;
}
#else
namespace lt {

template <typename C>
LT_NODISCARD constexpr auto data(C& c) -> decltype(c.data())
{
    return c.data();
}

template <typename C>
LT_NODISCARD constexpr auto data(C const& c) -> decltype(c.data())
{
    return c.data();
}

template <typename T, std::size_t N>
LT_NODISCARD constexpr auto data(T (&array)[N]) noexcept -> T*
{
    return array;
}

template <typename E>
LT_NODISCARD constexpr auto data(std::initializer_list<E> il) noexcept -> E const*
{
    return il.begin();
}
} // namespace lt
#endif