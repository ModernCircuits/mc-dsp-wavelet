#include "wavelib.h"

#include "helper.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <memory>

auto main() -> int
{
    auto obj = wavelet { "bior3.1" };

    auto const rows = 64;
    auto const cols = 48;
    auto const N = rows * cols;

    auto inp = makeZeros<double>(N);
    auto oup = makeZeros<double>(N);
    auto diff = makeZeros<double>(N);
    auto const J = 2;

    wt2_set* wt = wt2_init(obj, "swt", rows, cols, J);

    for (auto i = 0; i < rows; ++i) {
        for (auto k = 0; k < cols; ++k) {
            inp[i * cols + k] = generate_rnd();
            oup[i * cols + k] = 0.0;
        }
    }

    auto wavecoeffs = swt2(wt, inp.get());

    int ir { 0 };
    int ic { 0 };
    auto* cLL = getWT2Coeffs(wt, wavecoeffs.get(), J, "A", &ir, &ic);

    dispWT2Coeffs(cLL, ir, ic);

    iswt2(wt, wavecoeffs.get(), oup.get());

    for (auto i = 0; i < N; ++i) {
        diff[i] = oup[i] - inp[i];
    }

    wt2_summary(wt);
    std::printf("Abs Max %g \n", absmax(diff.get(), N));

    wt2_free(wt);

    return 0;
}
