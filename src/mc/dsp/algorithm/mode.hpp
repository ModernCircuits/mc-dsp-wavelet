// SPDX-License-Identifier: BSL-1.0

#pragma once

#include <mc/core/algorithm.hpp>
#include <mc/core/cassert.hpp>
#include <mc/core/cstddef.hpp>
#include <mc/core/span.hpp>

namespace mc::dsp {

template<typename T>
[[maybe_unused]] auto mode(Span<float const> v) -> float
{
    MC_ASSERT(ranges::is_sorted(v));

    auto const n   = v.size();
    float count    = 1;
    float countmax = 0;
    float current  = v[0];
    float moda     = 0;
    for (std::size_t i = 1; i < n; i++) {
        if (v[i] == current) {
            count++;
        } else if (count > countmax) {
            countmax = count;
            count    = 1;
            moda     = v[i - 1];
        }
        current = v[i];
    }
    return moda;
}

}  // namespace mc::dsp
