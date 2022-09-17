#pragma once

#include <mc/core/cmath.hpp>

namespace mc::dsp {

auto corrcoef(int n, float const* x, float const* y) -> float
{
    auto xm = 0.0F;
    auto ym = 0.0F;
    for (auto i = 0; i < n; ++i) {
        xm += x[i];
        ym += y[i];
    }

    xm = xm / n;
    ym = ym / n;

    auto num  = 0.0F;
    auto den1 = 0.0F;
    auto den2 = 0.0F;

    for (auto i = 0; i < n; ++i) {
        auto tx = x[i] - xm;
        auto ty = y[i] - ym;
        num += (tx * ty);
        den1 += (tx * tx);
        den2 += (ty * ty);
    }

    return num / sqrt(den1 * den2);
}

}  // namespace mc::dsp
