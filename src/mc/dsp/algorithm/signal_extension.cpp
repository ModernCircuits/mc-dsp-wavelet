// SPDX-License-Identifier: BSL-1.0

#include "signal_extension.hpp"

#include <mc/core/algorithm.hpp>

namespace mc {

auto toString(SignalExtension ext) -> String
{
    return ext == SignalExtension::periodic ? "periodic" : "symmetric";
}

template<typename T>
static auto
periodicExtensionImpl(T const* signal, std::size_t len, std::size_t a, T* output)
    -> std::size_t
{
    std::copy(signal, signal + len, output + a);

    auto len2 = len;
    if ((len % 2) != 0) {
        len2            = len + 1;
        output[a + len] = signal[len - 1];
    }

    for (std::size_t i = 0; i < a; ++i) {
        auto const temp1     = output[a + i];
        auto const temp2     = output[a + len2 - 1 - i];
        output[a - 1 - i]    = temp2;
        output[len2 + a + i] = temp1;
    }

    return len2;
}

template<typename T>
static auto
symmetricExtensionImpl(T const* signal, std::size_t len, std::size_t a, T* output)
    -> std::size_t
{
    // output is of length len + 2 * a
    std::copy(signal, signal + len, output + a);
    auto const newLength = len;
    for (std::size_t i = 0; i < a; ++i) {
        auto const temp1          = output[a + i];
        auto const temp2          = output[a + newLength - 1 - i];
        output[a - 1 - i]         = temp1;
        output[newLength + a + i] = temp2;
    }
    return newLength;
}

auto periodicExtension(Span<float const> in, std::size_t a, float* out) -> std::size_t
{
    return periodicExtensionImpl<float>(in.data(), in.size(), a, out);
}

auto periodicExtension(Span<double const> in, std::size_t a, double* out) -> std::size_t
{
    return periodicExtensionImpl<double>(in.data(), in.size(), a, out);
}

auto symmetricExtension(Span<float const> in, std::size_t a, float* out) -> std::size_t
{
    return symmetricExtensionImpl<float>(in.data(), in.size(), a, out);
}

auto symmetricExtension(Span<double const> in, std::size_t a, double* out) -> std::size_t
{
    return symmetricExtensionImpl<double>(in.data(), in.size(), a, out);
}

}  // namespace mc
