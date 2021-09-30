#include "lt/dsp/wavelets.hpp"

#include "lt/cmath.hpp"
#include "lt/cstdlib.hpp"
#include "lt/cstring.hpp"
#include "lt/format.hpp"
#include "lt/memory.hpp"

#include "lt/testing/test.hpp"
#include "readFileToVector.hpp"

namespace dsp = lt::dsp;

auto main() -> int
{
    auto obj = dsp::Wavelet { "bior3.5" };
    auto const n = std::size_t { 256 };

    auto const inp = readFileToVector("testData/signal.txt");

    auto out = std::make_unique<float[]>(n);
    auto diff = std::make_unique<float[]>(n);

    auto const j = std::size_t { 1 };

    auto wt = dsp::WaveletTransform(obj, "swt", n, j); // Initialize the wavelet transform object
    wt.convMethod(dsp::ConvolutionMethod::direct);

    swt(wt, inp.data()); // Perform SWT
    //SWT output can be accessed using wt.output vector. Use wt_summary to find out how to extract appx and detail coefficients

    for (std::size_t i = 0; i < wt.outlength; ++i) {
        fmt::printf("%g ", wt.output()[i]);
    }

    iswt(wt, out.get()); // Perform ISWT (if needed)
    // Test Reconstruction

    for (std::size_t i = 0; i < wt.signalLength(); ++i) {
        diff[i] = out[i] - inp[i];
    }

    fmt::printf("\n MAX %g \n", absmax(diff.get(), wt.signalLength())); // If Reconstruction succeeded then the output should be a small value.

    summary(wt); // Prints the full summary.

    return EXIT_SUCCESS;
}
