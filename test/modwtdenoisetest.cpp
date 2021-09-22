#include "wavelets/Denoise.hpp"

#include "readFileToVector.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>

static auto rmse(int n, double const* x, double const* y) -> double
{
    double rms;
    int i;

    rms = 0.0;

    for (i = 0; i < n; ++i) {
        rms += (x[i] - y[i]) * (x[i] - y[i]);
    }

    rms = std::sqrt(rms / (double)n);

    return rms;
}

static auto corrcoef(int n, double const* x, double const* y) -> double
{
    double cc;
    double xm;
    double ym;
    double tx;
    double ty;
    double num;
    double den1;
    double den2;
    int i;
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

    printf("Signal - Noisy Signal Stats \n");
    printf("RMSE %g\n", rmse(n, sig.data(), inp.data()));
    printf("Corr Coeff %g\n", corrcoef(n, sig.data(), inp.data()));

    printf("Signal - DeNoised Signal Stats \n");
    printf("RMSE %g\n", rmse(n, sig.data(), out.get()));
    printf("Corr Coeff %g\n", corrcoef(n, sig.data(), out.get()));

    return 0;
}
