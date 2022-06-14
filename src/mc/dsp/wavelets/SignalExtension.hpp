#pragma once

#include "mc/cstddef.hpp"
#include "mc/preprocessor.hpp"
#include "mc/string.hpp"

namespace mc::dsp
{

enum struct SignalExtension
{
    periodic,
    symmetric,
};

[[nodiscard]] inline auto toString(SignalExtension ext) -> std::string
{
    if (ext == SignalExtension::periodic) { return "periodic"; }
    return "symmetric";
}

template<typename T>
auto periodicExtension(T const* signal, std::size_t len, std::size_t a, T* output) -> std::size_t
{
    std::copy(signal, signal + len, output + a);

    auto len2 = len;
    if ((len % 2) != 0)
    {
        len2            = len + 1;
        output[a + len] = signal[len - 1];
    }

    for (std::size_t i = 0; i < a; ++i)
    {
        auto const temp1     = output[a + i];
        auto const temp2     = output[a + len2 - 1 - i];
        output[a - 1 - i]    = temp2;
        output[len2 + a + i] = temp1;
    }

    return len2;
}

template<typename T>
auto symmetricExtension(T const* signal, std::size_t len, std::size_t a, T* output) -> std::size_t
{
    // output is of length len + 2 * a
    std::copy(signal, signal + len, output + a);
    auto const newLength = len;
    for (std::size_t i = 0; i < a; ++i)
    {
        auto const temp1          = output[a + i];
        auto const temp2          = output[a + newLength - 1 - i];
        output[a - 1 - i]         = temp1;
        output[newLength + a + i] = temp2;
    }
    return newLength;
}
}  // namespace mc::dsp
