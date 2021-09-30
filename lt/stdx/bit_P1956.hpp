#pragma once

#if __has_include(<bit>)
#include <bit>
#endif

#include "lt/preprocessor.hpp"

#if defined(__cpp_lib_int_pow2)
namespace lt {
using std::bit_ceil;
using std::bit_floor;
using std::bit_width;
using std::has_single_bit;
}
#else

#include "lt/stdx/bit_P0553.hpp"

namespace lt {

template <typename T>
LT_NODISCARD constexpr auto bit_width(T x) noexcept // NOLINT(readability-identifier-naming)
    -> std::enable_if_t<detail::bitUnsignedInt<T>, int>
{
    return std::numeric_limits<T>::digits - countl_zero(x);
}

namespace detail {
    template <typename T, bool NoPromotion>
    struct BitCeilImpl;

    template <typename T>
    struct BitCeilImpl<T, false> {
        static auto ceil(T x) noexcept
        {
            // for types subject to integral promotion
            auto o = std::numeric_limits<unsigned>::digits - std::numeric_limits<T>::digits;
            return T { 1U << (bit_width(T { x - 1U }) + o) >> o };
        }
    };
    template <typename T>
    struct BitCeilImpl<T, true> {
        static auto ceil(T x) noexcept
        {
            return T { 1U } << bit_width(T { x - 1U });
        }
    };
}

template <typename T>
LT_NODISCARD constexpr auto bit_ceil(T x) noexcept // NOLINT(readability-identifier-naming)
    -> std::enable_if_t<detail::bitUnsignedInt<T>, T>
{
    if (x <= T(1)) {
        return T(1);
    }
    using impl_t = detail::BitCeilImpl<T, !std::is_same<T, decltype(+x)>::value>;
    return impl_t::ceil(x);
}
} // namespace lt
#endif
