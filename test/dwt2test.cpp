#include <mc/dsp/algorithm.hpp>
#include <mc/dsp/wavelet.hpp>

#include <mc/core/cmath.hpp>
#include <mc/core/cstdlib.hpp>
#include <mc/core/cstring.hpp>
#include <mc/core/format.hpp>
#include <mc/testing/test.hpp>

using namespace mc;

auto main() -> int
{
    auto obj = dsp::Wavelet{"db2"};

    auto rows = std::size_t{32};
    auto cols = std::size_t{30};
    auto n    = rows * cols;

    auto inp  = makeZeros<float>(n);
    auto oup  = makeZeros<float>(n);
    auto diff = makeZeros<float>(n);

    auto const j = std::size_t{3};

    auto wt = dsp::WaveletTransform2D(obj, "dwt", rows, cols, j);

    for (std::size_t i = 0; i < rows; ++i) {
        for (std::size_t k = 0; k < cols; ++k) {
            inp[i * cols + k] = generateRnd();
            oup[i * cols + k] = 0.0F;
        }
    }

    auto wavecoeffs = dwt(wt, inp.get());

    int ir{0};
    int ic{0};
    auto* cLL = dsp::getWT2Coeffs(wt, wavecoeffs.get(), 1, "D", &ir, &ic);

    dsp::dispWT2Coeffs(cLL, ir, ic);

    idwt(wt, wavecoeffs.get(), oup.get());

    for (std::size_t i = 0; i < rows * cols; ++i) { diff[i] = oup[i] - inp[i]; }

    print("{0}\n", wt);
    print("Abs Max {} \n", dsp::absmax(diff.get(), rows * cols));

    return EXIT_SUCCESS;
}
