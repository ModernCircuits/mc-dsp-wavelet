// SPDX-License-Identifier: BSL-1.0

#include "fft.hpp"

namespace mc {
auto makeFFT(std::size_t size) -> FFT<float>
{
    return FFT<float>{PFFFT_Complex_Float{size}};
}
}  // namespace mc
