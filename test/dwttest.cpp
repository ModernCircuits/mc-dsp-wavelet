#include "mc/dsp/wavelets.hpp"

#include "mc/cmath.hpp"
#include "mc/cstdlib.hpp"
#include "mc/cstring.hpp"
#include "mc/format.hpp"
#include "mc/memory.hpp"
#include "mc/utility.hpp"

#include "mc/testing/test.hpp"
#include "readFileToVector.hpp"

namespace dsp = mc::dsp;

auto main() -> int
{
    auto const input = readFileToVector("testData/signal.txt");
    auto const n     = 256;

    auto obj = dsp::Wavelet{"db4"};
    auto wt  = dsp::WaveletTransform(obj, "dwt", n, 3);
    wt.extension(dsp::SignalExtension::symmetric);
    wt.convMethod(dsp::ConvolutionMethod::direct);

    // DWT output can be accessed using wt.output vector.
    // Use wt_summary to find out how to extract appx and detail coefficients
    dwt(wt, input.data());

    for (auto i = 0; mc::cmp_less(i, wt.outlength); ++i) { fmt::printf("%g ", wt.output()[i]); }

    auto out = std::make_unique<float[]>(n);
    idwt(wt, out.get());

    auto diff = std::make_unique<float[]>(n);
    for (auto i = 0; mc::cmp_less(i, wt.signalLength()); ++i) { diff[i] = out[i] - input[i]; }

    // If Reconstruction succeeded then the output should be a small value.
    fmt::printf("\n MAX %g \n", absmax(diff.get(), wt.signalLength()));

    // Prints the full summary.
    summary(wt);

    return 0;
}
