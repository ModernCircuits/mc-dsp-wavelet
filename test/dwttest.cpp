#include "lt/dsp/wavelets.hpp"

#include "lt/cmath.hpp"
#include "lt/cstdlib.hpp"
#include "lt/cstring.hpp"
#include "lt/format.hpp"
#include "lt/memory.hpp"
#include "lt/utility.hpp"

#include "lt/testing/test.hpp"
#include "readFileToVector.hpp"

namespace dsp = lt::dsp;

auto main() -> int
{
    auto const input = readFileToVector("testData/signal.txt");
    auto const n = 256;

    auto obj = dsp::Wavelet { "db4" };
    auto wt = dsp::WaveletTransform(obj, "dwt", n, 3);
    wt.extension(dsp::SignalExtension::symmetric);
    wt.convMethod(dsp::ConvolutionMethod::direct);

    // DWT output can be accessed using wt.output vector.
    // Use wt_summary to find out how to extract appx and detail coefficients
    dwt(wt, input.data());

    for (auto i = 0; lt::cmp_less(i, wt.outlength); ++i) {
        fmt::printf("%g ", wt.output()[i]);
    }

    auto out = std::make_unique<float[]>(n);
    idwt(wt, out.get());

    auto diff = std::make_unique<float[]>(n);
    for (auto i = 0; lt::cmp_less(i, wt.signalLength()); ++i) {
        diff[i] = out[i] - input[i];
    }

    // If Reconstruction succeeded then the output should be a small value.
    fmt::printf("\n MAX %g \n", absmax(diff.get(), wt.signalLength()));

    // Prints the full summary.
    summary(wt);

    return 0;
}
