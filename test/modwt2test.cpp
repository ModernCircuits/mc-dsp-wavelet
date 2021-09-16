#include "wavelib.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>

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

auto generate_rnd() -> double
{
    double rnd;

    rnd = (double)(rand() % 100 + 1);

    return rnd;
}

auto main() -> int
{
    wave_object obj;
    wt2_object wt;
    int i;
    int k;
    int J;
    int rows;
    int cols;
    int N;
    int ir;
    int ic;
    double* inp;
    double* wavecoeffs;
    double* oup;
    double* cLL;
    double* diff;
    double amax;
    rows = 51;
    cols = 40;
    N = rows * cols;
    char const* name = "db2";
    obj = wave_init(name); // Initialize the wavelet

    inp = (double*)calloc(N, sizeof(double));
    oup = (double*)calloc(N, sizeof(double));
    diff = (double*)calloc(N, sizeof(double));
    J = 2;

    wt = wt2_init(obj, "modwt", rows, cols, J);

    for (i = 0; i < rows; ++i) {
        for (k = 0; k < cols; ++k) {
            //inp[i*cols + k] = i*cols + k;
            inp[i * cols + k] = generate_rnd();
            oup[i * cols + k] = 0.0;
        }
    }

    wavecoeffs = modwt2(wt, inp);

    cLL = getWT2Coeffs(wt, wavecoeffs, J, "A", &ir, &ic);

    //dispWT2Coeffs(cLL, ir, ic);

    imodwt2(wt, wavecoeffs, oup);

    for (i = 0; i < N; ++i) {
        diff[i] = oup[i] - inp[i];
    }

    amax = absmax(diff, N);

    wt2_summary(wt);

    printf("Abs Max %g \n", amax);

    wave_free(obj);
    wt2_free(wt);
    free(inp);
    free(wavecoeffs);
    free(oup);
    free(diff);
    return 0;
}
