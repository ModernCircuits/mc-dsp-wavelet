#include "lt/cmath.hpp"
#include "lt/dsp/wavelets.hpp"

#include "lt/testing/test.hpp"

#include <cstdio>
#include <cstdlib>
#include <cstring>

auto main() -> int
{
    auto obj = Wavelet { "db2" };

    auto rows = 32;
    auto cols = 30;
    auto n = rows * cols;

    auto inp = makeZeros<double>(n);
    auto oup = makeZeros<double>(n);
    auto diff = makeZeros<double>(n);

    auto const j = 3;

    auto wt = WaveletTransform2D(obj, "dwt", rows, cols, j);

    for (auto i = 0; i < rows; ++i) {
        for (auto k = 0; k < cols; ++k) {
            inp[i * cols + k] = generateRnd();
            oup[i * cols + k] = 0.0;
        }
    }

    auto wavecoeffs = dwt(wt, inp.get());

    int ir { 0 };
    int ic { 0 };
    auto* cLL = getWT2Coeffs(wt, wavecoeffs.get(), 1, "D", &ir, &ic);

    dispWT2Coeffs(cLL, ir, ic);

    idwt(wt, wavecoeffs.get(), oup.get());

    for (auto i = 0; i < rows * cols; ++i) {
        diff[i] = oup[i] - inp[i];
    }

    summary(wt);
    std::printf("Abs Max %g \n", absmax(diff.get(), rows * cols));

    return EXIT_SUCCESS;
}
