#pragma once

#include "lt/string.hpp"

namespace lt::dsp {

enum struct WindowFunction {
    rectangular,
    triangular,
    hann,
    hamming,
    blackman,
    blackmanHarris,
};

[[nodiscard]] auto toString(WindowFunction wf) -> std::string;

} // namespace lt::dsp
