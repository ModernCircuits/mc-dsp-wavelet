#include "wavelib.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>

auto absmax(double* array, int N) -> double
{
    double max;
    int i;

    max = 0.0;
    for (i = 0; i < N; ++i) {
        if (fabs(array[i]) >= max) {
            max = fabs(array[i]);
        }
    }

    return max;
}

auto main() -> int
{
    wave_object obj;
    wt_object wt;

    FILE* ifp;
    double temp[1200];

    char const* name = "db4";
    obj = wave_init(name); // Initialize the wavelet

    ifp = fopen("testData/signal.txt", "r");
    auto i = 0;
    if (ifp == nullptr) {
        printf("Cannot Open File");
        exit(100);
    }
    while (feof(ifp) == 0) {
        fscanf(ifp, "%lf \n", &temp[i]);
        i++;
    }
    auto N = 256;

    fclose(ifp);

    auto inp = std::make_unique<double[]>(N);
    auto out = std::make_unique<double[]>(N);
    auto diff = std::make_unique<double[]>(N);

    for (i = 0; i < N; ++i) {
        inp[i] = temp[i];
    }

    wt = wt_init(obj, "dwt", N, 3);
    setDWTExtension(wt, "sym");
    setWTConv(wt, "direct");

    // DWT output can be accessed using wt->output vector.
    // Use wt_summary to find out how to extract appx and detail coefficients
    dwt(wt, inp.get());

    for (i = 0; i < wt->outlength; ++i) {
        printf("%g ", wt->output[i]);
    }

    idwt(wt, out.get());
    for (i = 0; i < wt->siglength; ++i) {
        diff[i] = out[i] - inp[i];
    }

    // If Reconstruction succeeded then the output should be a small value.
    printf("\n MAX %g \n", absmax(diff.get(), wt->siglength));

    // Prints the full summary.
    wt_summary(wt);
    wave_free(obj);
    wt_free(wt);

    return 0;
}
