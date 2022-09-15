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

[[nodiscard]] auto toString(WindowFunction wf) -> std::string;
}  // namespace mc::dsp
