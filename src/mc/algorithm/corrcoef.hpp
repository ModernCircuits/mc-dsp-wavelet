// SPDX-License-Identifier: BSL-1.0

#pragma once

#include <mc/core/cstddef.hpp>

namespace mc {

[[nodiscard]] auto corrcoef(int n, float const* x, float const* y) -> float;

}  // namespace mc
