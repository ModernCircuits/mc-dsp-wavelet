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
    char const* name = "bior3.5";
    auto* obj = wave_init(name); // Initialize the wavelet
    auto const N = 256;

    auto const inp = readFileToVector("testData/signal.txt");

    auto out = std::make_unique<double[]>(N);
    auto diff = std::make_unique<double[]>(N);
    //wmean = mean(temp, N);

    auto const J = 1;

    wt_object wt = wt_init(obj, "swt", N, J); // Initialize the wavelet transform object
    setWTConv(wt, "direct");

    swt(wt, inp.data()); // Perform SWT
    //SWT output can be accessed using wt->output vector. Use wt_summary to find out how to extract appx and detail coefficients

    for (auto i = 0; i < wt->outlength; ++i) {
        std::printf("%g ", wt->output[i]);
    }

    iswt(wt, out.get()); // Perform ISWT (if needed)
    // Test Reconstruction

    for (auto i = 0; i < wt->siglength; ++i) {
        diff[i] = out[i] - inp[i];
    }

    std::printf("\n MAX %g \n", absmax(diff.get(), wt->siglength)); // If Reconstruction succeeded then the output should be a small value.

    wt_summary(wt); // Prints the full summary.

    wave_free(obj);
    wt_free(wt);

    return EXIT_SUCCESS;
}
