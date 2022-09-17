#pragma once

#include <mc/dsp/wavelet.hpp>

namespace mc::dsp {

struct TempoDetect
{
    TempoDetect(std::size_t n, std::size_t levels);

    [[nodiscard]] auto operator()(Span<float> input, float sampleRate) -> float;

private:
    dsp::Wavelet _wave;
    dsp::WaveletTransform _wt;
};

}  // namespace mc::dsp
