#include "lt/cmath.hpp"
#include "lt/dsp/wavelets.hpp"

#include "lt/testing/test.hpp"
#include "readFileToVector.hpp"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>

auto main() -> int
{
    auto const input = readFileToVector("testData/signal.txt");
    auto const n = 256;

    auto obj = Wavelet { "db4" };
    auto wt = WaveletTransform(obj, "dwt", n, 3);
    wt.extension(SignalExtension::symmetric);
    wt.convMethod(ConvolutionMethod::direct);

    // DWT output can be accessed using wt.output vector.
    // Use wt_summary to find out how to extract appx and detail coefficients
    dwt(wt, input.data());

    for (auto i = 0; i < wt.outlength; ++i) {
        printf("%g ", wt.output()[i]);
    }

    auto out = std::make_unique<double[]>(n);
    idwt(wt, out.get());

    auto diff = std::make_unique<double[]>(n);
    for (auto i = 0; i < wt.signalLength(); ++i) {
        diff[i] = out[i] - input[i];
    }

    // If Reconstruction succeeded then the output should be a small value.
    printf("\n MAX %g \n", absmax(diff.get(), wt.signalLength()));

    // Prints the full summary.
    summary(wt);

    return 0;
}
