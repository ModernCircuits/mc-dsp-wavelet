// SPDX-License-Identifier: BSL-1.0

#pragma once

#include <mc/dsp/convolution/convolute.hpp>
#include <mc/dsp/fft/rfft.hpp>

#include <mc/core/memory.hpp>
#include <mc/core/vector.hpp>

namespace mc::dsp {
struct FFTConvolver
{
    using value_type = float;

    FFTConvolver(std::size_t signalSize, std::size_t patchSize);

    auto convolute(Span<float const> signal, Span<float const> patch, float* output)
        -> void;

private:
    std::size_t _signalSize;
    std::size_t _patchSize;
    std::size_t _totalSize;

    RFFT<float> _fft;

    Vector<float> _signalScratch{};
    Vector<Complex<float>> _signalScratchOut{};

    Vector<float> _patchScratch{};
    Vector<Complex<float>> _patchScratchOut{};

    Vector<Complex<float>> _tmp{};
    Vector<float> _tmpOut{};
};
}  // namespace mc::dsp
