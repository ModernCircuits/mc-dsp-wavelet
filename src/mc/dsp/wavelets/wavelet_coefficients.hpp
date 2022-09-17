#pragma once

#include <mc/core/cstddef.hpp>
#include <mc/core/span.hpp>
#include <mc/core/string_view.hpp>

namespace mc::dsp {

template<typename T>
struct WaveletCoefficients
{
    StringView name;
    Span<T> coefficients;
    std::size_t length;
};

}  // namespace mc::dsp
