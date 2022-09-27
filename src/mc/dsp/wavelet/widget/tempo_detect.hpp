// SPDX-License-Identifier: BSL-1.0

#pragma once

#include <mc/dsp/wavelet.hpp>

namespace mc {

struct TempoDetect
{
    TempoDetect(std::size_t n, std::size_t levels);

    [[nodiscard]] auto operator()(Span<float> input, float sampleRate) -> float;

private:
    Wavelet _wave;
    WaveletTransform _wt;
};

}  // namespace mc
