// SPDX-License-Identifier: BSL-1.0

#pragma once

#include <mc/core/algorithm.hpp>
#include <mc/core/cstddef.hpp>

namespace mc {

template<typename T>
auto downSample(T const* x, std::size_t lenx, std::size_t m, T* y) -> void
{
    if (m == 0) { std::copy(x, x + static_cast<std::size_t>(lenx), y); }

    auto const n = (lenx - 1U) / m + 1U;
    for (std::size_t i = 0; i < n; ++i) { y[i] = x[i * m]; }
}

}  // namespace mc
