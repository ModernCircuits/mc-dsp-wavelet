#include "lt/cmath.hpp"
#include "lt/format.hpp"

#include "lt/dsp/wavelets.hpp"

#include "lt/testing/test.hpp"

#include "lt/cstdlib.hpp"
#include "lt/cstring.hpp"

auto main() -> int
{
    auto obj = Wavelet { "db2" };

    auto rows = std::size_t { 32 };
    auto cols = std::size_t { 30 };
    auto n = rows * cols;

    auto inp = makeZeros<float>(n);
    auto oup = makeZeros<float>(n);
    auto diff = makeZeros<float>(n);

    auto const j = std::size_t { 3 };

    auto wt = WaveletTransform2D(obj, "dwt", rows, cols, j);

    for (std::size_t i = 0; i < rows; ++i) {
        for (std::size_t k = 0; k < cols; ++k) {
            inp[i * cols + k] = generateRnd();
            oup[i * cols + k] = 0.0F;
        }
    }

    auto wavecoeffs = dwt(wt, inp.get());

    int ir { 0 };
    int ic { 0 };
    auto* cLL = getWT2Coeffs(wt, wavecoeffs.get(), 1, "D", &ir, &ic);

    dispWT2Coeffs(cLL, ir, ic);

    idwt(wt, wavecoeffs.get(), oup.get());

    for (std::size_t i = 0; i < rows * cols; ++i) {
        diff[i] = oup[i] - inp[i];
    }

    summary(wt);
    fmt::printf("Abs Max %g \n", absmax(diff.get(), rows * cols));

    return EXIT_SUCCESS;
}
