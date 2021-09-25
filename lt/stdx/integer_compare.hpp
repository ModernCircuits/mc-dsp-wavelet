#pragma once

#include <type_traits>
#include <utility>

#if not defined(__cpp_lib_integer_comparison_functions)
namespace std { // NOLINT(cert-dcl58-cpp)
template <typename T, typename U>
constexpr auto cmp_equal(T t, U u) noexcept -> bool // NOLINT(readability-identifier-naming)
{
    using UT = std::make_unsigned_t<T>;
    using UU = std::make_unsigned_t<U>;
    if constexpr (std::is_signed_v<T> == std::is_signed_v<U>) {
        return t == u;
    } else if constexpr (std::is_signed_v<T>) {
        return t < 0 ? false : UT(t) == u;
    } else {
        return u < 0 ? false : t == UU(u);
    }
}

template <typename T, typename U>
constexpr auto cmp_not_equal(T t, U u) noexcept -> bool // NOLINT(readability-identifier-naming)
{
    return !cmp_equal(t, u);
}

template <typename T, typename U>
constexpr auto cmp_less(T t, U u) noexcept -> bool // NOLINT(readability-identifier-naming)
{
    using UT = std::make_unsigned_t<T>;
    using UU = std::make_unsigned_t<U>;
    if constexpr (std::is_signed_v<T> == std::is_signed_v<U>) {
        return t < u;
    } else if constexpr (std::is_signed_v<T>) {
        return t < 0 ? true : UT(t) < u;
    } else {
        return u < 0 ? false : t < UU(u);
    }
}

template <typename T, typename U>
constexpr auto cmp_greater(T t, U u) noexcept -> bool // NOLINT(readability-identifier-naming)
{
    return cmp_less(u, t);
}

template <typename T, typename U>
constexpr auto cmp_less_equal(T t, U u) noexcept -> bool // NOLINT(readability-identifier-naming)
{
    return !cmp_greater(t, u);
}

template <typename T, typename U>
constexpr auto cmp_greater_equal(T t, U u) noexcept -> bool // NOLINT(readability-identifier-naming)
{
    return !cmp_less(t, u);
}
}
#endif