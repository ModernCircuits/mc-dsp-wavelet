#pragma once

#include "mc/dsp/wavelets.hpp"

namespace mc::dsp
{

struct TempoDetect
{
    TempoDetect(std::size_t n, std::size_t levels);

    [[nodiscard]] auto operator()(Span<float> input, float sampleRate) -> float;

private:
    dsp::Wavelet wave_;
    dsp::WaveletTransform wt_;
};

}  // namespace mc::dsp
