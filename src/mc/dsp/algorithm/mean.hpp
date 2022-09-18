#pragma once

#include <mc/core/numeric.hpp>
#include <mc/core/span.hpp>
#include <mc/core/utility.hpp>

namespace mc::dsp {

template<typename T>
auto mean(Span<T> range) -> T
{
    auto const sum = ranges::accumulate(range, T{0});
    auto const len = std::distance(ranges::begin(range), ranges::end(range));
    return sum / static_cast<T>(len);
}

}  // namespace mc::dsp
