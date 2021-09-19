#include "wauxlib.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>

static auto rmse(int N, double const* x, double const* y) -> double
{
    double rms;
    int i;

    rms = 0.0;

    for (i = 0; i < N; ++i) {
        rms += (x[i] - y[i]) * (x[i] - y[i]);
    }

    rms = sqrt(rms / (double)N);

    return rms;
}

static auto corrcoef(int N, double const* x, double const* y) -> double
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
    for (i = 0; i < N; ++i) {
        xm += x[i];
        ym += y[i];
    }

    xm = xm / N;
    ym = ym / N;
    num = den1 = den2 = 0.0;

    for (i = 0; i < N; ++i) {
        tx = x[i] - xm;
        ty = y[i] - ym;
        num += (tx * ty);
        den1 += (tx * tx);
        den2 += (ty * ty);
    }

    cc = num / sqrt(den1 * den2);

    return cc;
}

auto main() -> int
{
    // gcc -Wall -I../header -L../Bin modwtdenoisetest.c -o modwtdenoise -lwauxlib -lwavelib -lm
    /*
    modwtshrink can also be called from the denoise object. See denoisetest.c for more information
    */
    int i;
    int N;
    int J;
    FILE* ifp;

    double temp[2400];

    auto const* wname = "db5";
    auto const* ext = "per"; // The other option sym is only available with "fft" cmethod
    auto const* thresh = "soft";
    auto const* cmethod = "direct"; // The other option is "fft"

    ifp = fopen("testData/pieceregular1024.txt", "r");
    i = 0;
    if (ifp == nullptr) {
        printf("Cannot Open File");
        exit(100);
    }

    while (feof(ifp) == 0) {
        fscanf(ifp, "%lf \n", &temp[i]);
        i++;
    }

    fclose(ifp);

    N = i;
    J = 4;

    auto sig = std::make_unique<double[]>(N);
    auto inp = std::make_unique<double[]>(N);
    auto out = std::make_unique<double[]>(N);

    for (i = 0; i < N; ++i) {
        sig[i] = temp[i];
    }

    ifp = fopen("testData/PieceRegular10.txt", "r");
    i = 0;
    if (ifp == nullptr) {
        printf("Cannot Open File");
        exit(100);
    }

    while (feof(ifp) == 0) {
        fscanf(ifp, "%lf \n", &temp[i]);
        i++;
    }

    fclose(ifp);

    for (i = 0; i < N; ++i) {
        inp[i] = temp[i];
    }

    modwtshrink(sig.get(), N, J, wname, cmethod, ext, thresh, out.get());

    printf("Signal - Noisy Signal Stats \n");
    printf("RMSE %g\n", rmse(N, sig.get(), inp.get()));
    printf("Corr Coeff %g\n", corrcoef(N, sig.get(), inp.get()));

    printf("Signal - DeNoised Signal Stats \n");
    printf("RMSE %g\n", rmse(N, sig.get(), out.get()));
    printf("Corr Coeff %g\n", corrcoef(N, sig.get(), out.get()));

    return 0;
}
