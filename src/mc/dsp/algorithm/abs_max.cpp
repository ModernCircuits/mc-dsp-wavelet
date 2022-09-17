#pragma once

#include <mc/core/algorithm.hpp>
#include <mc/core/cassert.hpp>
#include <mc/core/cmath.hpp>
#include <mc/core/span.hpp>

namespace mc::dsp {

auto absmax(float const* array, std::size_t n) -> float
{
    MC_ASSERT(n != 0UL);

    auto const view = Span<float const>{array, n};
    auto const max  = ranges::max_element(view, {}, [](auto x) { return std::abs(x); });

    MC_ASSERT(max != ranges::end(view));
    return *max;
}
}  // namespace mc::dsp
