// SPDX-License-Identifier: BSL-1.0

#include <mc/dsp/fft.hpp>

#include <mc/core/array.hpp>
#include <mc/core/complex.hpp>
#include <mc/core/cstdlib.hpp>
#include <mc/core/iostream.hpp>

using namespace mc;

auto main() -> int
{
    auto engine = makeRFFT(4);
    auto in     = Array<float, 4>{0.0F, 1.0F, 0.3F, 0.0F};
    auto out    = Array<Complex<float>, 4>{};

    rfft(engine, in, out);

    for (auto c : out) { std::cout << c.real() << ' '; }
    std::cout << '\n';

    for (auto c : out) { std::cout << c.imag() << ' '; }
    std::cout << '\n';

    irfft(engine, out, in);
    for (auto c : in) { std::cout << c << ' '; }
    std::cout << '\n';

    return EXIT_SUCCESS;
}
