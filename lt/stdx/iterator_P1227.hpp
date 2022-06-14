#pragma once

#include <iterator>
#include <type_traits>

#include "lt/preprocessor.hpp"

#if defined(__cpp_lib_ssize)
namespace lt
{
using std::ssize;
}
#else
namespace lt
{
template<typename C>
LT_NODISCARD constexpr auto ssize(C const& c)
    -> std::common_type_t<std::ptrdiff_t, std::make_signed_t<decltype(c.size())> >
{
    using R = std::common_type_t<std::ptrdiff_t, std::make_signed_t<decltype(c.size())> >;
    return static_cast<R>(c.size());
}
template<typename T, std::ptrdiff_t N>
LT_NODISCARD constexpr auto ssize(T const (&array)[N]) noexcept -> std::ptrdiff_t
{
    (void)array;
    return N;
}
}  // namespace lt
#endif