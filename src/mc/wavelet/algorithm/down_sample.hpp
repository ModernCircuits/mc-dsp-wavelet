// SPDX-License-Identifier: BSL-1.0

#pragma once

#include <mc/core/algorithm.hpp>
#include <mc/core/cstddef.hpp>

namespace mc {

template<typename T>
auto downSample(T const* x, size_t lenx, size_t m, T* y) -> void
{
    if (m == 0) { std::copy(x, x + static_cast<size_t>(lenx), y); }

    auto const n = (lenx - 1U) / m + 1U;
    for (size_t i = 0; i < n; ++i) { y[i] = x[i * m]; }
}

}  // namespace mc
