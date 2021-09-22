#include "lt/dsp/wavelets.hpp"

#include "helper.hpp"
#include "readFileToVector.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>

auto main() -> int
{
    auto obj = Wavelet { "db4" };
    summary(obj);

    auto n = 177;
    auto j = 2;

    auto inp = readFileToVector("testData/signal.txt");
    auto out = std::make_unique<double[]>(n);
    auto diff = std::make_unique<double[]>(n);

    // Initialize the wavelet transform object
    auto wt = WaveletTransform(obj, "modwt", n, j);

    // MODWT output can be accessed using wt.output vector.
    // Use wt_summary to find out how to extract appx and detail coefficients
    modwt(wt, inp.data()); // Perform MODWT

    for (auto i = 0; i < wt.outlength; ++i) {
        printf("%g ", wt.output()[i]);
    }

    imodwt(wt, out.get()); // Perform ISWT (if needed)

    for (auto i = 0; i < wt.siglength; ++i) {
        diff[i] = out[i] - inp[i];
    }

    // If Reconstruction succeeded then the output should be a small value.
    printf("\n MAX %g \n", absmax(diff.get(), wt.siglength));

    summary(wt); // Prints the full summary.

    return 0;
}
