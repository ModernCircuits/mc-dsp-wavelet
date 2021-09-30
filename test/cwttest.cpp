#include "lt/dsp/wavelets.hpp"

#include "lt/cmath.hpp"
#include "lt/format.hpp"
#include "lt/memory.hpp"

#include "readFileToVector.hpp"

namespace dsp = lt::dsp;

auto main() -> int
{
    int i = 0;
    int n = 0;
    int j = 0;
    int subscale = 0;
    int a0 = 0;
    int iter = 0;
    int nd = 0;
    int k = 0;
    float dt = NAN;
    float dj = NAN;
    float s0 = NAN;
    float param = NAN;
    float mn = NAN;
    float td = NAN;
    float tn = NAN;
    float den = NAN;
    float num = NAN;
    float reconMean = NAN;
    float reconVar = NAN;

    // Set Morlet wavelet. Other options "paul" and "dog"
    auto const* wave = "morlet";
    auto const* type = "pow";

    n = 504;
    param = 6.0F;
    subscale = 4;
    dt = 0.25;
    s0 = dt;
    dj = 1.0F / (float)subscale;
    j = 11 * subscale; // Total Number of scales
    a0 = 2; //power

    auto wt = dsp::ContinuousWaveletTransform { wave, param, n, dt, j };

    auto const inp = readFileToVector("testData/sst_nino3.dat");
    auto oup = std::make_unique<float[]>(n);

    wt.scales(s0, dj, type, a0);

    cwt(wt, inp.data());

    fmt::printf("\n MEAN %g \n", wt.smean);

    mn = 0.0F;

    for (i = 0; i < n; ++i) {
        mn += std::sqrt(wt.output[i].real() * wt.output[i].real() + wt.output[i].imag() * wt.output[i].imag());
    }

    summary(wt);

    fmt::printf("\n abs mean %g \n", mn / n);

    fmt::printf("\n\n");
    fmt::printf("Let CWT w = w(j, n/2 - 1) where n = %d\n\n", n);
    nd = n / 2 - 1;

    fmt::printf("%-15s%-15s%-15s%-15s \n", "j", "Scale", "Period", "ABS(w)^2");
    for (k = 0; k < wt.J; ++k) {
        iter = nd + k * n;
        fmt::printf("%-15d%-15lf%-15lf%-15lf \n", k, wt.scale[k], wt.period[k],
            wt.output[iter].real() * wt.output[iter].real() + wt.output[iter].imag() * wt.output[iter].imag());
    }

    icwt(wt, oup.get());

    num = den = reconMean = 0.0F;
    fmt::printf("\n\n");
    fmt::printf("Signal Reconstruction\n");
    fmt::printf("%-15s%-15s%-15s \n", "i", "Input(i)", "Output(i)");

    for (i = n - 10; i < n; ++i) {
        fmt::printf("%-15d%-15lf%-15lf \n", i, inp[i], oup[i]);
    }

    for (i = 0; i < n; ++i) {
        //fmt::printf("%g %g \n", oup[i] ,inp[i] - wt.smean);
        td = inp[i];
        tn = oup[i] - td;
        num += (tn * tn);
        den += (td * td);
        reconMean += oup[i];
    }

    reconVar = std::sqrt(num / n);
    reconMean /= n;

    fmt::printf("\nRMS Error %g \n", std::sqrt(num) / std::sqrt(den));
    fmt::printf("\nVariance %g \n", reconVar);
    fmt::printf("\nMean %g \n", reconMean);

    return 0;
}
