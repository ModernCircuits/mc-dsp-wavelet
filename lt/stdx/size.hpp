
#pragma once

#include <iterator>
#include <type_traits>

#if defined(__cpp_lib_nonmember_container_access)
namespace lt {
using std::size;
}
#else
namespace lt {
template <typename C>
[[nodiscard]] constexpr auto size(const C& c) -> decltype(c.size())
{
    return c.size();
}
template <typename T, std::size_t N>
[[nodiscard]] constexpr auto size(const T (&array)[N]) noexcept -> std::size_t
{
    return N;
}
} // namespace lt
#endif