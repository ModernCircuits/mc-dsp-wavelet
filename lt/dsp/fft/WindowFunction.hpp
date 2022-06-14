#pragma once

#include "lt/preprocessor.hpp"
#include "lt/string.hpp"

namespace lt
{
namespace dsp
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
}  // namespace dsp
}  // namespace lt
