// SPDX-License-Identifier: BSL-1.0

#include "WindowFunction.hpp"

#include <mc/core/cassert.hpp>

namespace mc::dsp {
auto toString(WindowFunction wf) -> String
{
    switch (wf) {
        case WindowFunction::rectangular: return "Rectangular";
        case WindowFunction::triangular: return "Triangular";
        case WindowFunction::hann: return "Hann";
        case WindowFunction::hamming: return "Hamming";
        case WindowFunction::blackman: return "Blackman";
        case WindowFunction::blackmanHarris: return "Blackman-Harris";
    }

    MC_ASSERT(false);  // NOLINT
    return "";
}
}  // namespace mc::dsp
