// SPDX-License-Identifier: BSL-1.0

#include <mc/dsp/wavelet.hpp>

#include <mc/core/cstdlib.hpp>
#include <mc/core/print.hpp>

using namespace mc;

auto main() -> int
{
    print("{0}\n", Wavelet{"db1"});
    print("{0}\n", Wavelet{"db2"});
    print("{0}\n", Wavelet{"db3"});
    return EXIT_SUCCESS;
}
