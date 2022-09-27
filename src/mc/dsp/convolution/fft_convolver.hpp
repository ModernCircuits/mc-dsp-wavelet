// SPDX-License-Identifier: BSL-1.0

#pragma once

#include <mc/dsp/convolution/convolute.hpp>
#include <mc/dsp/fft/rfft.hpp>

#include <mc/core/memory.hpp>
#include <mc/core/vector.hpp>

namespace mc {
struct FFTConvolver
{
    using value_type = float;

    FFTConvolver(size_t signalSize, size_t patchSize);

    auto convolute(Span<float const> signal, Span<float const> patch, float* output)
        -> void;

private:
    size_t _signalSize;
    size_t _patchSize;
    size_t _totalSize;

    RFFT<float> _fft;

    Vector<float> _signalScratch{};
    Vector<Complex<float>> _signalScratchOut{};

    Vector<float> _patchScratch{};
    Vector<Complex<float>> _patchScratchOut{};

    Vector<Complex<float>> _tmp{};
    Vector<float> _tmpOut{};
};
}  // namespace mc
