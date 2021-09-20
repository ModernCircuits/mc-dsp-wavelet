#include "wavelib.h"

#include "helper.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>

auto main() -> int
{
    auto obj = wavelet { "db4" };
    auto const N = 788 + 23;
    auto const J = 4;

    auto inp = std::make_unique<double[]>(N);
    auto oup = std::make_unique<double[]>(N);
    auto diff = std::make_unique<double[]>(N);

    for (auto i = 1; i < N + 1; ++i) {
        inp[i - 1] = i;
    }

    wpt_set* wt = wpt_init(&obj, N, J);
    setDWPTExtension(wt, "per");
    setDWPTEntropy(wt, "logenergy", 0);

    dwpt(wt, inp.get());

    idwpt(wt, oup.get());

    for (auto i = 0; i < N; ++i) {
        diff[i] = (inp[i] - oup[i]) / inp[i];
    }

    wpt_summary(wt);

    // If Reconstruction succeeded then the output should be a small value.
    printf("\n MAX %g \n", absmax(diff.get(), wt->siglength));

    wpt_free(wt);
    return 0;
}
