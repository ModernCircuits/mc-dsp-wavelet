#include <mc/dsp/wavelets.hpp>

#include <mc/core/cmath.hpp>
#include <mc/core/cstdlib.hpp>
#include <mc/core/cstring.hpp>
#include <mc/core/format.hpp>
#include <mc/core/memory.hpp>
#include <mc/core/utility.hpp>
#include <mc/testing/test.hpp>

using namespace mc;

auto main() -> int
{
    auto const input = readFileToVector<float>("test_data/raw/signal.txt");
    auto const n     = 256;

    auto obj = dsp::Wavelet{"db4"};
    auto wt  = dsp::WaveletTransform(obj, "dwt", n, 3);
    wt.extension(dsp::SignalExtension::symmetric);
    wt.convMethod(dsp::ConvolutionMethod::direct);

    // DWT output can be accessed using wt.output vector.
    // Use wt_summary to find out how to extract appx and detail coefficients
    dwt(wt, input.data());

    for (auto i = 0; mc::cmp_less(i, wt.outlength); ++i) { print("{} ", wt.output()[i]); }

    auto out = makeUnique<float[]>(n);
    idwt(wt, out.get());

    auto diff = makeUnique<float[]>(n);
    for (auto i = 0; mc::cmp_less(i, wt.signalLength()); ++i) {
        diff[i] = out[i] - input[i];
    }

    // If Reconstruction succeeded then the output should be a small value.
    print("\n MAX {} \n", absmax(diff.get(), wt.signalLength()));

    // Prints the full summary.
    summary(wt);

    return 0;
}
