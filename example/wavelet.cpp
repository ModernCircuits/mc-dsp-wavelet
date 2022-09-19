// SPDX-License-Identifier: BSL-1.0

#include <mc/dsp/wavelet.hpp>

#include <mc/core/cstdlib.hpp>
#include <mc/core/format.hpp>

using namespace mc;

auto main() -> int
{
    print("{0}\n", dsp::Wavelet{"db1"});
    print("{0}\n", dsp::Wavelet{"db2"});
    print("{0}\n", dsp::Wavelet{"db3"});
    return EXIT_SUCCESS;
}
