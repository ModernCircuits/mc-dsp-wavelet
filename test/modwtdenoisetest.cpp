#include "mc/dsp/wavelets/Denoise.hpp"

#include "mc/cmath.hpp"
#include "mc/cstdlib.hpp"
#include "mc/cstring.hpp"
#include "mc/format.hpp"
#include "mc/memory.hpp"

#include "mc/testing/test.hpp"

namespace dsp = mc::dsp;

static auto rmse(int n, float const* x, float const* y) -> float
{
    auto rms = 0.0F;
    for (std::size_t i = 0; i < static_cast<std::size_t>(n); ++i) { rms += (x[i] - y[i]) * (x[i] - y[i]); }
    return std::sqrt(rms / (float)n);
}

static auto corrcoef(int n, float const* x, float const* y) -> float
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
    for (i = 0; i < n; ++i)
    {
        xm += x[i];
        ym += y[i];
    }

    xm  = xm / n;
    ym  = ym / n;
    num = den1 = den2 = 0.0F;

    for (i = 0; i < n; ++i)
    {
        tx = x[i] - xm;
        ty = y[i] - ym;
        num += (tx * ty);
        den1 += (tx * tx);
        den2 += (ty * ty);
    }

    cc = num / std::sqrt(den1 * den2);

    return cc;
}

// modwtshrink can also be called from the denoise object.
// See denoisetest.cpp for more information.
auto main() -> int
{
    auto const* wname   = "db5";
    auto const* ext     = "per";
    auto const* thresh  = "soft";
    auto const* cmethod = "direct";

    auto const inp = readFileToVector<float>("testData/PieceRegular10.txt");
    auto sig       = readFileToVector<float>("testData/pieceregular1024.txt");

    auto const n = sig.size();
    auto const j = 4;
    auto out     = std::make_unique<float[]>(n);

    dsp::modwtshrink(sig.data(), n, j, wname, cmethod, ext, thresh, out.get());

    fmt::printf("Signal - Noisy Signal Stats \n");
    fmt::printf("RMSE %g\n", rmse(n, sig.data(), inp.data()));
    fmt::printf("Corr Coeff %g\n", corrcoef(n, sig.data(), inp.data()));

    fmt::printf("Signal - DeNoised Signal Stats \n");
    fmt::printf("RMSE %g\n", rmse(n, sig.data(), out.get()));
    fmt::printf("Corr Coeff %g\n", corrcoef(n, sig.data(), out.get()));

    return 0;
}
