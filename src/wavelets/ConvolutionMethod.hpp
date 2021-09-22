#pragma once

#include <string>

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