#include "wavelets.hpp"

#include "helper.hpp"
#include "readFileToVector.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>

auto main() -> int
{
    auto obj = Wavelet { "bior3.5" };
    auto const n = 256;

    auto const inp = readFileToVector("testData/signal.txt");

    auto out = std::make_unique<double[]>(n);
    auto diff = std::make_unique<double[]>(n);

    auto const j = 1;

    auto wt = WaveletTransform(obj, "swt", n, j); // Initialize the wavelet transform object
    wt.convMethod(ConvolutionMethod::direct);

    swt(wt, inp.data()); // Perform SWT
    //SWT output can be accessed using wt.output vector. Use wt_summary to find out how to extract appx and detail coefficients

    for (auto i = 0; i < wt.outlength; ++i) {
        std::printf("%g ", wt.output()[i]);
    }

    iswt(wt, out.get()); // Perform ISWT (if needed)
    // Test Reconstruction

    for (auto i = 0; i < wt.siglength; ++i) {
        diff[i] = out[i] - inp[i];
    }

    std::printf("\n MAX %g \n", absmax(diff.get(), wt.siglength)); // If Reconstruction succeeded then the output should be a small value.

    summary(wt); // Prints the full summary.

    return EXIT_SUCCESS;
}
