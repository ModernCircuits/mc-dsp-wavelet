#include "lt/dsp/wavelets.hpp"

#include "lt/cmath.hpp"
#include "lt/cstdlib.hpp"
#include "lt/cstring.hpp"
#include "lt/format.hpp"
#include "lt/memory.hpp"

#include "lt/testing/test.hpp"

auto main() -> int
{
    auto obj = Wavelet { "bior3.1" };

    auto const rows = std::size_t { 64 };
    auto const cols = std::size_t { 48 };
    auto const n = rows * cols;

    auto inp = makeZeros<double>(n);
    auto oup = makeZeros<double>(n);
    auto diff = makeZeros<double>(n);
    auto const j = std::size_t { 2 };

    auto wt = WaveletTransform2D(obj, "swt", rows, cols, j);

    for (std::size_t i = 0; i < rows; ++i) {
        for (std::size_t k = 0; k < cols; ++k) {
            inp[i * cols + k] = generateRnd();
            oup[i * cols + k] = 0.0;
        }
    }

    auto wavecoeffs = swt2(wt, inp.get());

    int ir { 0 };
    int ic { 0 };
    auto* cLL = getWT2Coeffs(wt, wavecoeffs.get(), j, "A", &ir, &ic);

    dispWT2Coeffs(cLL, ir, ic);

    iswt2(wt, wavecoeffs.get(), oup.get());

    for (std::size_t i = 0; i < n; ++i) {
        diff[i] = oup[i] - inp[i];
    }

    summary(wt);
    fmt::printf("Abs Max %g \n", absmax(diff.get(), n));

    return 0;
}
