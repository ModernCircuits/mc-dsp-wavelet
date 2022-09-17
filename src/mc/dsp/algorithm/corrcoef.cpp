#pragma once

#include <mc/core/cmath.hpp>

namespace mc::dsp {

auto corrcoef(int n, float const* x, float const* y) -> float
{
    float cc   = NAN;
    float xm   = NAN;
    float ym   = NAN;
    float tx   = NAN;
    float ty   = NAN;
    float num  = NAN;
    float den1 = NAN;
    float den2 = NAN;
    int i      = 0;
    xm = ym = 0.0F;
    for (i = 0; i < n; ++i) {
        xm += x[i];
        ym += y[i];
    }

    xm  = xm / n;
    ym  = ym / n;
    num = den1 = den2 = 0.0F;

    for (i = 0; i < n; ++i) {
        tx = x[i] - xm;
        ty = y[i] - ym;
        num += (tx * ty);
        den1 += (tx * tx);
        den2 += (ty * ty);
    }

    cc = num / sqrt(den1 * den2);

    return cc;
}

}  // namespace mc::dsp
