#include "wavelets.hpp"

#include <cmath>
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

    WaveletTree* wt = wtreeInit(&obj, n, j);
    setWTREEExtension(wt, "sym");

    wtree(wt, inp.get());
    summary(*wt);
    auto const x = 3;
    auto const y = 5;
    auto const len = getWTREENodelength(wt, x);
    std::printf("\n %d", len);
    std::printf("\n");
    auto oup = std::make_unique<double[]>(len);

    std::printf("Node [%d %d] Coefficients : \n", x, y);
    getWTREECoeffs(wt, x, y, oup.get(), len);
    for (auto i = 0; i < len; ++i) {
        std::printf("%g ", oup[i]);
    }
    std::printf("\n");

    wtreeFree(wt);
    return EXIT_SUCCESS;
}
