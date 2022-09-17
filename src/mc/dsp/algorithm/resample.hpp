#pragma once

#include <mc/core/algorithm.hpp>
#include <mc/core/cstddef.hpp>
#include <mc/core/utility.hpp>

namespace mc::dsp {

template<typename T>
auto upSample(T const* x, std::size_t lenx, std::size_t m, T* y) -> void
{
    if (m == 0) { std::copy(x, x + static_cast<std::size_t>(lenx), y); }

    auto j       = std::size_t{1};
    auto k       = std::size_t{0};
    auto const n = m * (lenx - 1U) + 1U;
    for (std::size_t i = 0; i < n; ++i) {
        j--;
        y[i] = T(0);
        if (j == 0) {
            y[i] = x[k];
            k++;
            j = m;
        }
    }
}

// Returns even numbered output. Last value is set to zero
template<typename SrcIt, typename DestIt>
auto upSampleEven(SrcIt srcF, SrcIt srcL, DestIt destF, std::size_t m) -> void
{
    using T = typename std::iterator_traits<SrcIt>::value_type;
    if (m == 0) { std::copy(srcF, srcL, destF); }

    auto it   = destF;
    auto last = destF + (m * static_cast<std::size_t>(std::distance(srcF, srcL)));
    auto j    = std::size_t{1};
    for (; it != last; ++it) {
        j--;
        *it = T(0);
        if (j == 0) {
            *it = *srcF;
            ++srcF;
            j = m;
        }
    }
}

template<typename T>
auto downSample(T const* x, std::size_t lenx, std::size_t m, T* y) -> void
{
    if (m == 0) { std::copy(x, x + static_cast<std::size_t>(lenx), y); }

    auto const n = (lenx - 1U) / m + 1U;
    for (std::size_t i = 0; i < n; ++i) { y[i] = x[i * m]; }
}

}  // namespace mc::dsp
