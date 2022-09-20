// SPDX-License-Identifier: BSL-1.0

#include "rfft.hpp"

namespace mc::dsp {
auto makeRFFT(std::size_t size) -> RFFT<float>
{
    return RFFT<float>{PFFFT_Real_Float{size}};
}
}  // namespace mc::dsp
