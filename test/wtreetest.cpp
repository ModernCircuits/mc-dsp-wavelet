#include "wavelib.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>

auto main() -> int
{
    auto obj = wavelet { "db3" };
    auto const N = 147;
    auto inp = std::make_unique<double[]>(N);
    for (auto i = 1; i < N + 1; ++i) {
        inp[i - 1] = -0.25 * i * i * i + 25 * i * i + 10 * i;
    }
    auto const J = 3;

    wtree_set* wt = wtree_init(&obj, N, J);
    setWTREEExtension(wt, "sym");

    wtree(wt, inp.get());
    wtree_summary(wt);
    auto const X = 3;
    auto const Y = 5;
    auto const len = getWTREENodelength(wt, X);
    std::printf("\n %d", len);
    std::printf("\n");
    auto oup = std::make_unique<double[]>(len);

    std::printf("Node [%d %d] Coefficients : \n", X, Y);
    getWTREECoeffs(wt, X, Y, oup.get(), len);
    for (auto i = 0; i < len; ++i) {
        std::printf("%g ", oup[i]);
    }
    std::printf("\n");

    wtree_free(wt);
    return EXIT_SUCCESS;
}
