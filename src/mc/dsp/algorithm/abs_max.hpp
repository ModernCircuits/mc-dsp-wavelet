#pragma once

#include <mc/core/cstddef.hpp>

namespace mc::dsp {

[[nodiscard]] auto absmax(float const* array, std::size_t n) -> float;

}  // namespace mc::dsp
