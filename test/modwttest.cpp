#include "wavelib.h"

#include "helper.hpp"
#include "readFileToVector.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>

auto main() -> int
{
    auto obj = wavelet { "db4" };
    wave_summary(obj);

    auto N = 177;
    auto J = 2;

    auto inp = readFileToVector("testData/signal.txt");
    auto out = std::make_unique<double[]>(N);
    auto diff = std::make_unique<double[]>(N);

    // Initialize the wavelet transform object
    auto wt = wavelet_transform(obj, "modwt", N, J);

    // MODWT output can be accessed using wt.output vector.
    // Use wt_summary to find out how to extract appx and detail coefficients
    modwt(&wt, inp.data()); // Perform MODWT

    for (auto i = 0; i < wt.outlength; ++i) {
        printf("%g ", wt.output[i]);
    }

    imodwt(&wt, out.get()); // Perform ISWT (if needed)

    for (auto i = 0; i < wt.siglength; ++i) {
        diff[i] = out[i] - inp[i];
    }

    // If Reconstruction succeeded then the output should be a small value.
    printf("\n MAX %g \n", absmax(diff.get(), wt.siglength));

    wt_summary(&wt); // Prints the full summary.

    return 0;
}
