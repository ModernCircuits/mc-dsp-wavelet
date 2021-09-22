#include "wavelib.h"

#include "readFileToVector.hpp"

#include <cmath>
#include <memory>

auto main() -> int
{
    int i;
    int n;
    int j;
    int subscale;
    int a0;
    int iter;
    int nd;
    int k;
    double dt;
    double dj;
    double s0;
    double param;
    double mn;
    double td;
    double tn;
    double den;
    double num;
    double reconMean;
    double reconVar;
    ComplexWaveletTransform* wt;

    // Set Morlet wavelet. Other options "paul" and "dog"
    auto const* wave = "morlet";
    auto const* type = "pow";

    n = 504;
    param = 6.0;
    subscale = 4;
    dt = 0.25;
    s0 = dt;
    dj = 1.0 / (double)subscale;
    j = 11 * subscale; // Total Number of scales
    a0 = 2; //power

    wt = cwtInit(wave, param, n, dt, j);

    auto const inp = readFileToVector("testData/sst_nino3.dat");
    auto oup = std::make_unique<double[]>(n);

    setCWTScales(wt, s0, dj, type, a0);

    cwt(wt, inp.data());

    std::printf("\n MEAN %g \n", wt->smean);

    mn = 0.0;

    for (i = 0; i < n; ++i) {
        mn += std::sqrt(wt->output[i].re * wt->output[i].re + wt->output[i].im * wt->output[i].im);
    }

    summary(*wt);

    std::printf("\n abs mean %g \n", mn / n);

    std::printf("\n\n");
    std::printf("Let CWT w = w(j, n/2 - 1) where n = %d\n\n", n);
    nd = n / 2 - 1;

    std::printf("%-15s%-15s%-15s%-15s \n", "j", "Scale", "Period", "ABS(w)^2");
    for (k = 0; k < wt->J; ++k) {
        iter = nd + k * n;
        std::printf("%-15d%-15lf%-15lf%-15lf \n", k, wt->scale[k], wt->period[k],
            wt->output[iter].re * wt->output[iter].re + wt->output[iter].im * wt->output[iter].im);
    }

    icwt(wt, oup.get());

    num = den = reconMean = 0.0;
    std::printf("\n\n");
    std::printf("Signal Reconstruction\n");
    std::printf("%-15s%-15s%-15s \n", "i", "Input(i)", "Output(i)");

    for (i = n - 10; i < n; ++i) {
        std::printf("%-15d%-15lf%-15lf \n", i, inp[i], oup[i]);
    }

    for (i = 0; i < n; ++i) {
        //std::printf("%g %g \n", oup[i] ,inp[i] - wt->smean);
        td = inp[i];
        tn = oup[i] - td;
        num += (tn * tn);
        den += (td * td);
        reconMean += oup[i];
    }

    reconVar = std::sqrt(num / n);
    reconMean /= n;

    std::printf("\nRMS Error %g \n", std::sqrt(num) / std::sqrt(den));
    std::printf("\nVariance %g \n", reconVar);
    std::printf("\nMean %g \n", reconMean);

    cwtFree(wt);
    return 0;
}
