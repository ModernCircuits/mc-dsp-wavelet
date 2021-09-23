#include "lt/cmath.hpp"
#include "lt/dsp/wavelets.hpp"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>

auto main() -> int
{
    auto obj = Wavelet { "db3" };
    auto const n = 147;
    auto inp = std::make_unique<double[]>(n);
    for (auto i = 1; i < n + 1; ++i) {
        inp[i - 1] = -0.25 * i * i * i + 25 * i * i + 10 * i;
    }
    auto const j = 3;

    auto wt = WaveletTree(&obj, n, j);
    wt.extension("sym");

    wtree(wt, inp.get());
    summary(wt);
    auto const x = 3;
    auto const y = 5;
    auto const len = wt.nodeLength(x);
    std::printf("\n %d", len);
    std::printf("\n");
    auto oup = std::make_unique<double[]>(len);

    std::printf("Node [%d %d] Coefficients : \n", x, y);
    wt.coeffs(x, y, oup.get(), len);
    for (auto i = 0; i < len; ++i) {
        std::printf("%g ", oup[i]);
    }
    std::printf("\n");

    return EXIT_SUCCESS;
}
