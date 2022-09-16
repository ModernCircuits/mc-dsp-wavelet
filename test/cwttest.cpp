#include <mc/dsp/wavelets.hpp>

#include <mc/core/cmath.hpp>
#include <mc/core/format.hpp>
#include <mc/core/memory.hpp>
#include <mc/testing/test.hpp>

#include <fmt/printf.h>

using namespace mc;

auto main() -> int
{
    // Set Morlet wavelet. Other options "paul" and "dog"
    auto const* wave = "morlet";
    auto const* type = "pow";

    auto const n        = 504;
    auto const param    = 6.0F;
    auto const subscale = 4;
    auto const dt       = 0.25;
    auto const s0       = dt;
    auto const dj       = 1.0F / (float)subscale;
    auto const j        = 11 * subscale;  // Total Number of scales
    auto const a0       = 2;              // power

    auto wt = dsp::ContinuousWaveletTransform{wave, param, n, dt, j};

    auto const inp = readFileToVector<float>("test_data/raw/sst_nino3.dat");
    auto oup       = makeUnique<float[]>(n);

    wt.scales(s0, dj, type, a0);

    cwt(wt, inp.data());

    print("\n MEAN {} \n", wt.smean);

    auto mn = 0.0F;

    for (auto i = 0; i < n; ++i) {
        mn += mc::sqrt(
            wt.output[i].real() * wt.output[i].real()
            + wt.output[i].imag() * wt.output[i].imag()
        );
    }

    summary(wt);

    print("\n abs mean {} \n", mn / n);

    print("\n\n");
    print("Let CWT w = w(j, n/2 - 1) where n = {}\n\n", n);
    auto nd = n / 2 - 1;

    fmt::printf("%-15s%-15s%-15s%-15s \n", "j", "Scale", "Period", "ABS(w)^2");
    for (auto k = 0; k < wt.J; ++k) {
        auto iter = nd + k * n;
        fmt::printf(
            "%-15d%-15lf%-15lf%-15lf \n",
            k,
            wt.scale[k],
            wt.period[k],
            wt.output[iter].real() * wt.output[iter].real()
                + wt.output[iter].imag() * wt.output[iter].imag()
        );
    }

    icwt(wt, oup.get());

    auto num       = 0.0F;
    auto den       = 0.0F;
    auto reconMean = 0.0F;
    print("\n\n");
    print("Signal Reconstruction\n");
    fmt::printf("%-15s%-15s%-15s \n", "i", "Input(i)", "Output(i)");

    for (auto i = n - 10; i < n; ++i) {
        fmt::printf("%-15d%-15lf%-15lf \n", i, inp[i], oup[i]);
    }

    for (auto i = 0; i < n; ++i) {
        auto const td = inp[i];
        auto const tn = oup[i] - td;
        num += (tn * tn);
        den += (td * td);
        reconMean += oup[i];
    }

    auto const reconVar = mc::sqrt(num / n);
    reconMean /= n;

    print("\nRMS Error {} \n", mc::sqrt(num) / mc::sqrt(den));
    print("\nVariance {} \n", reconVar);
    print("\nMean {} \n", reconMean);

    return 0;
}
