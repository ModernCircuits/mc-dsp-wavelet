// SPDX-License-Identifier: BSL-1.0

#include <mc/wavelet/algorithm.hpp>

#include <mc/core/cmath.hpp>
#include <mc/core/cstdlib.hpp>
#include <mc/core/cstring.hpp>
#include <mc/core/memory.hpp>
#include <mc/core/print.hpp>
#include <mc/testing/test.hpp>
#include <mc/wavelet.hpp>

using namespace mc;

auto main() -> int
{
    auto obj = Wavelet{"db4"};

    auto const rows = size_t{64};
    auto const cols = size_t{48};
    auto const n    = rows * cols;

    auto inp     = makeZeros<float>(n);
    auto oup     = makeZeros<float>(n);
    auto diff    = makeZeros<float>(n);
    auto const j = size_t{2};

    auto wt = WaveletTransform2D(obj, "swt", rows, cols, j);

    for (size_t i = 0; i < rows; ++i) {
        for (size_t k = 0; k < cols; ++k) {
            inp[i * cols + k] = generateRnd();
            oup[i * cols + k] = 0.0;
        }
    }

    auto wavecoeffs = swt2(wt, inp.get());

    int ir{0};
    int ic{0};
    auto* cLL = getWT2Coeffs(wt, wavecoeffs.get(), j, "A", &ir, &ic);

    dispWT2Coeffs(cLL, ir, ic);

    iswt2(wt, wavecoeffs.get(), oup.get());

    for (size_t i = 0; i < n; ++i) { diff[i] = oup[i] - inp[i]; }

    print("{0}\n", wt);
    print("Abs Max {} \n", absmax(diff.get(), n));

    return 0;
}
