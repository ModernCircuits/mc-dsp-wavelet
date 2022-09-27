// SPDX-License-Identifier: BSL-1.0

#pragma once

#include <mc/core/cstddef.hpp>

namespace mc {

[[nodiscard]] auto relError(float const* data, float const* rec, size_t n) -> float;

}  // namespace mc
