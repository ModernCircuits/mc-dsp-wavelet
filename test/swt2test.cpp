#include "lt/cmath.hpp"
#include "lt/dsp/wavelets.hpp"

#include "lt/testing/test.hpp"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <memory>

auto main() -> int
{
    auto obj = Wavelet { "bior3.1" };

    auto const rows = 64;
    auto const cols = 48;
    auto const n = rows * cols;

    auto inp = makeZeros<double>(n);
    auto oup = makeZeros<double>(n);
    auto diff = makeZeros<double>(n);
    auto const j = 2;

    auto wt = WaveletTransform2D(obj, "swt", rows, cols, j);

    for (auto i = 0; i < rows; ++i) {
        for (auto k = 0; k < cols; ++k) {
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

    for (auto i = 0; i < n; ++i) {
        diff[i] = oup[i] - inp[i];
    }

    summary(wt);
    std::printf("Abs Max %g \n", absmax(diff.get(), n));

    return 0;
}
