// SPDX-License-Identifier: BSL-1.0

#pragma once

#include <mc/wavelet.hpp>

namespace mc {

struct TempoDetect
{
    TempoDetect(size_t n, size_t levels);

    [[nodiscard]] auto operator()(Span<float> input, float sampleRate) -> float;

private:
    Wavelet _wave;
    WaveletTransform _wt;
};

}  // namespace mc
