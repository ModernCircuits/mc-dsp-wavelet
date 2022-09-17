#include <mc/dsp/wavelets.hpp>

#include <mc/core/cstdlib.hpp>

using namespace mc;

auto main() -> int
{
    dsp::summary(dsp::Wavelet{"db1"});
    dsp::summary(dsp::Wavelet{"db2"});
    dsp::summary(dsp::Wavelet{"db3"});
    return EXIT_SUCCESS;
}
