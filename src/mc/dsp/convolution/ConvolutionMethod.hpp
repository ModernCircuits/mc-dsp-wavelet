#pragma once

#include <mc/core/config.hpp>

#include <mc/core/string.hpp>

namespace mc::dsp {
enum struct ConvolutionMethod
{
    direct,
    fft,
};

[[nodiscard]] inline auto toString(ConvolutionMethod method) -> String
{
    if (method == ConvolutionMethod::direct) { return "direct"; }
    return "fft";
}
}  // namespace mc::dsp
