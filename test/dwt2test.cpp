#include "wavelib.h"

#include "helper.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

auto main() -> int
{
    wave_object obj = wave_init("db2");

    auto rows = 32;
    auto cols = 30;
    auto N = rows * cols;

    auto inp = makeZeros<double>(N);
    auto oup = makeZeros<double>(N);
    auto diff = makeZeros<double>(N);

    auto const J = 3;

    wt2_object wt = wt2_init(obj, "dwt", rows, cols, J);

    for (auto i = 0; i < rows; ++i) {
        for (auto k = 0; k < cols; ++k) {
            inp[i * cols + k] = generate_rnd();
            oup[i * cols + k] = 0.0;
        }
    }

    auto* wavecoeffs = dwt2(wt, inp.get());

    int ir { 0 };
    int ic { 0 };
    auto* cLL = getWT2Coeffs(wt, wavecoeffs, 1, "D", &ir, &ic);

    dispWT2Coeffs(cLL, ir, ic);

    idwt2(wt, wavecoeffs, oup.get());

    for (auto i = 0; i < rows * cols; ++i) {
        diff[i] = oup[i] - inp[i];
    }

    wt2_summary(wt);
    std::printf("Abs Max %g \n", absmax(diff.get(), rows * cols));

    wave_free(obj);
    wt2_free(wt);
    std::free(wavecoeffs);

    return EXIT_SUCCESS;
}
