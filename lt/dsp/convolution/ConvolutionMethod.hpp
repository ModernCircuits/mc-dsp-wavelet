#pragma once

#include "lt/string.hpp"

enum struct ConvolutionMethod {
    direct,
    fft,
};

[[nodiscard]] inline auto toString(ConvolutionMethod method) -> std::string
{
    if (method == ConvolutionMethod::direct) {
        return "direct";
    }
    return "fft";
}