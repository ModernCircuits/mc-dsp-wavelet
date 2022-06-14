#include "mc/dsp/wavelets.hpp"

#include "mc/cmath.hpp"
#include "mc/cstdlib.hpp"
#include "mc/cstring.hpp"
#include "mc/format.hpp"
#include "mc/memory.hpp"

namespace dsp = mc::dsp;

auto main() -> int
{
    auto obj     = dsp::Wavelet{"db3"};
    auto const n = 147;
    auto inp     = std::make_unique<float[]>(n);
    for (auto i = 1; i < n + 1; ++i) { inp[i - 1] = -0.25 * i * i * i + 25 * i * i + 10 * i; }
    auto const j = 3;

    auto wt = dsp::WaveletTree(&obj, n, j);
    wt.extension("sym");

    wtree(wt, inp.get());
    summary(wt);
    auto const x   = 3;
    auto const y   = 5;
    auto const len = wt.nodeLength(x);
    fmt::printf("\n %d", len);
    fmt::printf("\n");
    auto oup = std::make_unique<float[]>(len);

    fmt::printf("Node [%d %d] Coefficients : \n", x, y);
    wt.coeffs(x, y, oup.get(), len);
    for (std::size_t i = 0; i < len; ++i) { fmt::printf("%g ", oup[i]); }
    fmt::printf("\n");

    return EXIT_SUCCESS;
}
