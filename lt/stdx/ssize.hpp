#pragma once

#include <iterator>
#include <type_traits>

#if not defined(__cpp_lib_ssize)
namespace std { // NOLINT(cert-dcl58-cpp)
template <typename C>
[[nodiscard]] constexpr auto ssize(const C& c)
    -> std::common_type_t<std::ptrdiff_t,
        std::make_signed_t<decltype(c.size())>>
{
    using R = std::common_type_t<std::ptrdiff_t,
        std::make_signed_t<decltype(c.size())>>;
    return static_cast<R>(c.size());
}
template <typename T, std::ptrdiff_t N>
[[nodiscard]] constexpr auto ssize(const T (&array)[N]) noexcept -> std::ptrdiff_t
{
    (void)array;
    return N;
}
} // namespace std
#endif