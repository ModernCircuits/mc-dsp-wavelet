// SPDX-License-Identifier: BSL-1.0

#pragma once

#include <mc/core/config.hpp>

#include <mc/core/string.hpp>

namespace mc::dsp {
enum struct WindowFunction
{
    rectangular,
    triangular,
    hann,
    hamming,
    blackman,
    blackmanHarris,
};

[[nodiscard]] auto toString(WindowFunction wf) -> String;
}  // namespace mc::dsp
