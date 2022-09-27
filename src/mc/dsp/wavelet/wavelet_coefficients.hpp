// SPDX-License-Identifier: BSL-1.0

#pragma once

#include <mc/core/cstddef.hpp>
#include <mc/core/span.hpp>
#include <mc/core/string_view.hpp>

namespace mc {

template<typename T>
struct WaveletCoefficients
{
    StringView name;
    Span<T const> coefficients;
    size_t length{};
};

}  // namespace mc
