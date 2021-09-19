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
    wave_object obj = wave_init("db4");

    auto const input = readFileToVector("testData/signal.txt");
    auto const N = 256;

    wt_object wt = wt_init(obj, "dwt", N, 3);
    setDWTExtension(wt, "sym");
    setWTConv(wt, "direct");

    // DWT output can be accessed using wt->output vector.
    // Use wt_summary to find out how to extract appx and detail coefficients
    dwt(wt, input.data());

    for (auto i = 0; i < wt->outlength; ++i) {
        printf("%g ", wt->output[i]);
    }

    auto out = std::make_unique<double[]>(N);
    idwt(wt, out.get());

    auto diff = std::make_unique<double[]>(N);
    for (auto i = 0; i < wt->siglength; ++i) {
        diff[i] = out[i] - input[i];
    }

    // If Reconstruction succeeded then the output should be a small value.
    printf("\n MAX %g \n", absmax(diff.get(), wt->siglength));

    // Prints the full summary.
    wt_summary(wt);
    wave_free(obj);
    wt_free(wt);

    return 0;
}
