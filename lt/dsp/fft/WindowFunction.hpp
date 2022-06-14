#pragma once

#include "lt/preprocessor.hpp"
#include "lt/string.hpp"

namespace lt::dsp
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

LT_NODISCARD auto toString(WindowFunction wf) -> std::string;
}  // namespace lt::dsp
