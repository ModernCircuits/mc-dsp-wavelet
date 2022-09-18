#pragma once

#include <mc/core/algorithm.hpp>
#include <mc/core/span.hpp>
#include <mc/core/utility.hpp>

namespace mc::dsp {

template<typename T>
[[nodiscard]] auto indexOfPeak(Span<T> range) -> std::size_t
{
    auto [min, max] = ranges::minmax_element(range);
    if (std::abs(*min) >= std::abs(*max)) {
        return std::distance(ranges::begin(range), min);
    }
    return std::distance(ranges::begin(range), max);
}

}  // namespace mc::dsp
