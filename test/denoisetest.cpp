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

    int i;
    int N;
    int J;
    FILE* ifp;

    denoise_object obj;
    double temp[2400];

    char const* wname = "db5";
    char const* method = "dwt"; // Available - dwt, swt and modwt. modwt works only with modwtshrink. The other two methods work with
    // visushrink and sureshrink
    char const* ext = "sym"; // sym and per work with dwt. swt and modwt only use per extension when called through denoise.
    // You can use sy extension if you directly call modwtshrink with cmethod set to fft. See modwtdenoisetest.c file
    char const* thresh = "soft"; // soft or hard
    char const* level = "all"; // noise estimation at "first" or "all" levels. modwt only has the option of "all"

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

    auto inp = std::make_unique<double[]>(N);
    auto oup = std::make_unique<double[]>(N);
    auto sig = std::make_unique<double[]>(N);

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
    obj = denoise_init(N, J, wname);
    setDenoiseMethod(obj, "visushrink"); // sureshrink is also the default. The other option with dwt and swt is visushrink.
    // modwt works only with modwtshrink method
    setDenoiseWTMethod(obj, method); // Default is dwt. the other options are swt and modwt
    setDenoiseWTExtension(obj, ext); // Default is sym. the other option is per
    setDenoiseParameters(obj, thresh, level); // Default for thresh is soft. Other option is hard
    // Default for level is all. The other option is first

    denoise(obj, inp.get(), oup.get());

    // Alternative to denoise_object
    // Just use visushrink, modwtshrink and sureshrink functions
    //visushrink(inp.get(),N,J,wname,method,ext,thresh,level,oup.get());
    //sureshrink(inp.get(),N,J,wname,method,ext,thresh,level,oup.get());
    // modwtshrink(sig.get(),N,J,wname,cmethod,ext,thresh,oup.get()); See modwtdenoisetest.c
    //ofp = fopen("testData/denoiseds.txt", "w");

    printf("Signal - Noisy Signal Stats \n");
    printf("RMSE %g\n", rmse(N, sig.get(), inp.get()));
    printf("Corr Coeff %g\n", corrcoef(N, sig.get(), inp.get()));

    printf("Signal - DeNoised Signal Stats \n");
    printf("RMSE %g\n", rmse(N, sig.get(), oup.get()));
    printf("Corr Coeff %g\n", corrcoef(N, sig.get(), oup.get()));

    denoise_free(obj);
    return 0;
}
