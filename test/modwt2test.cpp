#include "wavelib.h"

#include "helper.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>

auto main() -> int
{
    wave_object obj = wave_init("db2"); // Initialize the wavelet

    auto rows = 51;
    auto cols = 40;
    auto N = rows * cols;

    auto inp = makeZeros<double>(N);
    auto oup = makeZeros<double>(N);
    auto diff = makeZeros<double>(N);

    auto J = 2;
    wt2_object wt = wt2_init(obj, "modwt", rows, cols, J);

    for (auto i = 0; i < rows; ++i) {
        for (auto k = 0; k < cols; ++k) {
            //inp[i*cols + k] = i*cols + k;
            inp[i * cols + k] = generate_rnd();
            oup[i * cols + k] = 0.0;
        }
    }

    auto* wavecoeffs = modwt2(wt, inp.get());

    int ir { 0 };
    int ic { 0 };
    auto* cLL = getWT2Coeffs(wt, wavecoeffs, J, "A", &ir, &ic);

    imodwt2(wt, wavecoeffs, oup.get());

    for (auto i = 0; i < N; ++i) {
        diff[i] = oup[i] - inp[i];
    }

    wt2_summary(wt);
    printf("Abs Max %g \n", absmax(diff.get(), N));

    wave_free(obj);
    wt2_free(wt);
    free(wavecoeffs);
    return 0;
}
