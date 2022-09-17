#pragma once

#include <mc/core/cstddef.hpp>

namespace mc::dsp {

[[nodiscard]] auto rmsError(float const* data, float const* rec, std::size_t n) -> float;

}  // namespace mc::dsp
