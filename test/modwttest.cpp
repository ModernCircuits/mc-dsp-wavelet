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
    for (auto i = 0; i < N; ++i) {
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
    obj = wave_init(name);
    wave_summary(obj);

    ifp = fopen("testData/signal.txt", "r");
    auto idx = 0;
    if (ifp == nullptr) {
        printf("Cannot Open File");
        exit(100);
    }
    while (feof(ifp) == 0) {
        fscanf(ifp, "%lf \n", &temp[idx]);
        ++idx;
    }
    auto N = 177;

    fclose(ifp);

    auto inp = std::make_unique<double[]>(N);
    auto out = std::make_unique<double[]>(N);
    auto diff = std::make_unique<double[]>(N);

    for (auto i = 0; i < N; ++i) {
        inp[i] = temp[i];
        //printf("%g \n",inp[i]);
    }
    auto J = 2;

    // Initialize the wavelet transform object
    wt = wt_init(obj, "modwt", N, J);

    modwt(wt, inp.get()); // Perform MODWT
    //MODWT output can be accessed using wt->output vector. Use wt_summary to find out how to extract appx and detail coefficients

    for (auto i = 0; i < wt->outlength; ++i) {
        printf("%g ", wt->output[i]);
    }

    imodwt(wt, out.get()); // Perform ISWT (if needed)
    // Test Reconstruction

    for (auto i = 0; i < wt->siglength; ++i) {
        diff[i] = out[i] - inp[i];
    }

    printf("\n MAX %g \n", absmax(diff.get(), wt->siglength)); // If Reconstruction succeeded then the output should be a small value.

    wt_summary(wt); // Prints the full summary.

    wave_free(obj);
    wt_free(wt);

    return 0;
}
