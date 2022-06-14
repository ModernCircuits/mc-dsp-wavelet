#pragma once

#include "mc/preprocessor.hpp"
#include "mc/string.hpp"

namespace mc::dsp
{
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
