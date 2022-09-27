// SPDX-License-Identifier: BSL-1.0

#pragma once

#include <mc/core/config.hpp>

#include <mc/core/string.hpp>

namespace mc {
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
}  // namespace mc
