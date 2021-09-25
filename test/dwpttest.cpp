#include "lt/dsp/wavelets.hpp"

#include "lt/cmath.hpp"
#include "lt/format.hpp"

#include "lt/testing/test.hpp"

#include "lt/cstdlib.hpp"
#include "lt/cstring.hpp"
#include "lt/memory.hpp"

auto main() -> int
{
    auto obj = Wavelet { "db4" };
    auto const n = std::size_t { 788 + 23 };
    auto const j = std::size_t { 4 };

    auto inp = std::make_unique<double[]>(n);
    auto oup = std::make_unique<double[]>(n);
    auto diff = std::make_unique<double[]>(n);

    for (std::size_t i = 1; i < n + 1; ++i) {
        inp[i - 1] = i;
    }

    auto wt = WaveletPacketTransform(&obj, n, j);
    setDWPTExtension(wt, "per");
    setDWPTEntropy(wt, "logenergy", 0);

    dwt(wt, inp.get());

    idwt(wt, oup.get());

    for (std::size_t i = 0; i < n; ++i) {
        diff[i] = (inp[i] - oup[i]) / inp[i];
    }

    summary(wt);

    // If Reconstruction succeeded then the output should be a small value.
    fmt::printf("\n MAX %g \n", absmax(diff.get(), wt.siglength));

    return 0;
}
