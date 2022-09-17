#include <mc/dsp/algorithm.hpp>
#include <mc/dsp/wavelet.hpp>

#include <mc/core/cmath.hpp>
#include <mc/core/cstdlib.hpp>
#include <mc/core/cstring.hpp>
#include <mc/core/format.hpp>
#include <mc/core/memory.hpp>
#include <mc/testing/test.hpp>

using namespace mc;

auto main() -> int
{
    auto obj = dsp::Wavelet{"db4"};

    auto const rows = std::size_t{64};
    auto const cols = std::size_t{48};
    auto const n    = rows * cols;

    auto inp     = makeZeros<float>(n);
    auto oup     = makeZeros<float>(n);
    auto diff    = makeZeros<float>(n);
    auto const j = std::size_t{2};

    auto wt = dsp::WaveletTransform2D(obj, "swt", rows, cols, j);

    for (std::size_t i = 0; i < rows; ++i) {
        for (std::size_t k = 0; k < cols; ++k) {
            inp[i * cols + k] = generateRnd();
            oup[i * cols + k] = 0.0;
        }
    }

    auto wavecoeffs = swt2(wt, inp.get());

    int ir{0};
    int ic{0};
    auto* cLL = getWT2Coeffs(wt, wavecoeffs.get(), j, "A", &ir, &ic);

    dsp::dispWT2Coeffs(cLL, ir, ic);

    iswt2(wt, wavecoeffs.get(), oup.get());

    for (std::size_t i = 0; i < n; ++i) { diff[i] = oup[i] - inp[i]; }

    print("{0}\n", summary(wt));
    print("Abs Max {} \n", dsp::absmax(diff.get(), n));

    return 0;
}
