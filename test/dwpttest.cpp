#include <mc/dsp/wavelets.hpp>

#include <mc/core/cmath.hpp>
#include <mc/core/cstdlib.hpp>
#include <mc/core/cstring.hpp>
#include <mc/core/format.hpp>
#include <mc/core/memory.hpp>
#include <mc/testing/test.hpp>

#include <fmt/printf.h>

using namespace mc;

auto main() -> int
{
    auto obj     = dsp::Wavelet{"db4"};
    auto const n = std::size_t{788 + 23};
    auto const j = std::size_t{4};

    auto inp  = makeUnique<float[]>(n);
    auto oup  = makeUnique<float[]>(n);
    auto diff = makeUnique<float[]>(n);

    for (std::size_t i = 1; i < n + 1; ++i) { inp[i - 1] = i; }

    auto wt = dsp::WaveletPacketTransform(&obj, n, j);
    setDWPTExtension(wt, "per");
    setDWPTEntropy(wt, "logenergy", 0);

    dwpt(wt, inp.get());

    idwpt(wt, oup.get());

    for (std::size_t i = 0; i < n; ++i) { diff[i] = (inp[i] - oup[i]) / inp[i]; }

    summary(wt);

    // If Reconstruction succeeded then the output should be a small value.
    fmt::printf("\n MAX %g \n", absmax(diff.get(), wt.signalLength()));

    return 0;
}
