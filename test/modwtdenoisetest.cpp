#include "lt/dsp/wavelets/Denoise.hpp"

#include "lt/cmath.hpp"
#include "lt/format.hpp"

#include "readFileToVector.hpp"

#include "lt/cstdlib.hpp"
#include "lt/cstring.hpp"
#include "lt/memory.hpp"

static auto rmse(int n, double const* x, double const* y) -> double
{
    double rms = NAN;
    int i = 0;

    rms = 0.0;

    for (i = 0; i < n; ++i) {
        rms += (x[i] - y[i]) * (x[i] - y[i]);
    }

    rms = std::sqrt(rms / (double)n);

    return rms;
}

static auto corrcoef(int n, double const* x, double const* y) -> double
{
    double cc = NAN;
    double xm = NAN;
    double ym = NAN;
    double tx = NAN;
    double ty = NAN;
    double num = NAN;
    double den1 = NAN;
    double den2 = NAN;
    int i = 0;
    xm = ym = 0.0;
    for (i = 0; i < n; ++i) {
        xm += x[i];
        ym += y[i];
    }

    xm = xm / n;
    ym = ym / n;
    num = den1 = den2 = 0.0;

    for (i = 0; i < n; ++i) {
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
    auto const* wname = "db5";
    auto const* ext = "per";
    auto const* thresh = "soft";
    auto const* cmethod = "direct";

    auto const inp = readFileToVector("testData/PieceRegular10.txt");
    auto sig = readFileToVector("testData/pieceregular1024.txt");

    auto const n = sig.size();
    auto const j = 4;
    auto out = std::make_unique<double[]>(n);

    modwtshrink(sig.data(), n, j, wname, cmethod, ext, thresh, out.get());

    fmt::printf("Signal - Noisy Signal Stats \n");
    fmt::printf("RMSE %g\n", rmse(n, sig.data(), inp.data()));
    fmt::printf("Corr Coeff %g\n", corrcoef(n, sig.data(), inp.data()));

    fmt::printf("Signal - DeNoised Signal Stats \n");
    fmt::printf("RMSE %g\n", rmse(n, sig.data(), out.get()));
    fmt::printf("Corr Coeff %g\n", corrcoef(n, sig.data(), out.get()));

    return 0;
}
