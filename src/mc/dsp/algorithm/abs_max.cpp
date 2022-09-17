#pragma once

#include <mc/core/algorithm.hpp>
#include <mc/core/cmath.hpp>

namespace mc::dsp {

auto absmax(float const* array, std::size_t n) -> float
{
    auto max = 0.0F;
    for (auto i = std::size_t{0}; i < n; ++i) {
        if (std::abs(array[i]) >= max) { max = std::abs(array[i]); }
    }
    return max;
}
}  // namespace mc::dsp
