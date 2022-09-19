// SPDX-License-Identifier: BSL-1.0

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

}  // namespace mc::dsp
