#pragma once

#include <mc/core/cmath.hpp>

namespace mc::dsp {

auto rmsError(float const* data, float const* rec, std::size_t n) -> float
{
    float sum = 0;
    for (std::size_t i = 0; i < n; ++i) { sum += (data[i] - rec[i]) * (data[i] - rec[i]); }
    return sqrt(sum / ((float)n - 1));
}

}  // namespace mc::dsp
