#pragma once

#include <type_traits>
#include <utility>

#if defined(__cpp_lib_integer_comparison_functions)
namespace lt {
using std::cmp_equal;
using std::cmp_greater;
using std::cmp_greater_equal;
using std::cmp_less;
using std::cmp_less_equal;
using std::cmp_not_equal;
} // namespace lt
#else
namespace lt {
template <typename T, typename U>
constexpr bool cmp_equal(T t, U u) noexcept
{
    using UT = std::make_unsigned_t<T>;
    using UU = std::make_unsigned_t<U>;
    if constexpr (std::is_signed_v<T> == std::is_signed_v<U>)
        return t == u;
    else if constexpr (std::is_signed_v<T>)
        return t < 0 ? false : UT(t) == u;
    else
        return u < 0 ? false : t == UU(u);
}

template <typename T, typename U>
constexpr bool cmp_not_equal(T t, U u) noexcept
{
    return !cmp_equal(t, u);
}

template <typename T, typename U>
constexpr bool cmp_less(T t, U u) noexcept
{
    using UT = std::make_unsigned_t<T>;
    using UU = std::make_unsigned_t<U>;
    if constexpr (std::is_signed_v<T> == std::is_signed_v<U>)
        return t < u;
    else if constexpr (std::is_signed_v<T>)
        return t < 0 ? true : UT(t) < u;
    else
        return u < 0 ? false : t < UU(u);
}

template <typename T, typename U>
constexpr bool cmp_greater(T t, U u) noexcept
{
    return cmp_less(u, t);
}

template <typename T, typename U>
constexpr bool cmp_less_equal(T t, U u) noexcept
{
    return !cmp_greater(t, u);
}

template <typename T, typename U>
constexpr bool cmp_greater_equal(T t, U u) noexcept
{
    return !cmp_less(t, u);
}
}
#endif