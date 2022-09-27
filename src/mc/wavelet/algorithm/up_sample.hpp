// SPDX-License-Identifier: BSL-1.0

#pragma once

#include <mc/core/algorithm.hpp>
#include <mc/core/cstddef.hpp>
#include <mc/core/utility.hpp>

namespace mc {

template<typename T>
auto upSample(T const* x, size_t lenx, size_t m, T* y) -> void
{
    if (m == 0) { std::copy(x, x + static_cast<size_t>(lenx), y); }

    auto j       = size_t{1};
    auto k       = size_t{0};
    auto const n = m * (lenx - 1U) + 1U;
    for (size_t i = 0; i < n; ++i) {
        j--;
        y[i] = T(0);
        if (j == 0) {
            y[i] = x[k];
            k++;
            j = m;
        }
    }
}

}  // namespace mc
