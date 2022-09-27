// SPDX-License-Identifier: BSL-1.0

#include <mc/dsp/wavelet.hpp>

#include <mc/core/cmath.hpp>
#include <mc/core/cstdlib.hpp>
#include <mc/core/cstring.hpp>
#include <mc/core/memory.hpp>
#include <mc/core/print.hpp>

using namespace mc;

auto main() -> int
{
    auto obj     = dsp::Wavelet{"db3"};
    auto const n = 147;
    auto inp     = makeUnique<float[]>(n);
    for (auto i = 1; i < n + 1; ++i) {
        inp[i - 1] = -0.25F * i * i * i + 25.0F * i * i + 10.0F * i;
    }
    auto const j = 3;

    auto wt = dsp::WaveletTree(&obj, n, j);
    wt.extension("sym");

    wtree(wt, inp.get());
    print("{0}\n", wt);
    auto const x   = 3;
    auto const y   = 5;
    auto const len = wt.nodeLength(x);
    print("\n {}\n", len);
    auto oup = makeUnique<float[]>(len);

    print("Node [{} {}] Coefficients : \n", x, y);
    wt.coeffs(x, y, oup.get(), len);
    for (std::size_t i = 0; i < len; ++i) { print("{} ", oup[i]); }
    print("\n");

    return EXIT_SUCCESS;
}
