// SPDX-License-Identifier: BSL-1.0

#pragma once

#include <mc/core/config.hpp>

#include <mc/core/cstddef.hpp>
#include <mc/core/span.hpp>
#include <mc/core/string.hpp>

namespace mc {

enum struct SignalExtension
{
    periodic,
    symmetric,
};

[[nodiscard]] auto toString(SignalExtension ext) -> String;

auto periodicExtension(Span<float const> in, size_t a, float* out) -> size_t;
auto periodicExtension(Span<double const> in, size_t a, double* out) -> size_t;

auto symmetricExtension(Span<float const> in, size_t a, float* out) -> size_t;
auto symmetricExtension(Span<double const> in, size_t a, double* out) -> size_t;

}  // namespace mc
