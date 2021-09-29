#pragma once

#include "lt/preprocessor.hpp"
#include "lt/string.hpp"

enum struct ConvolutionMethod {
    direct,
    fft,
};

LT_NODISCARD inline auto toString(ConvolutionMethod method) -> std::string
{
    if (method == ConvolutionMethod::direct) {
        return "direct";
    }
    return "fft";
}