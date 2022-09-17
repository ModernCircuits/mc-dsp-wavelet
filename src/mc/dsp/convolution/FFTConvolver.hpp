#pragma once

#include <mc/dsp/convolution/convolute.hpp>
#include <mc/dsp/fft/FFT.hpp>

namespace mc::dsp {
struct FFTConvolver
{
    using value_type = float;

    FFTConvolver(std::size_t signalSize, std::size_t patchSize);

    auto convolute(float const* signal, float const* patch, float* output) const -> void;

private:
    std::size_t _signalSize;
    std::size_t _patchSize;
    std::size_t _totalSize;

    UniquePtr<RFFT> _forwardFFT;
    UniquePtr<RFFT> _backwardFFT;

    UniquePtr<float[]> _signalScratch{makeUnique<float[]>(_totalSize)};
    UniquePtr<Complex<float>[]> _signalScratchOut {
        makeUnique<Complex<float>[]>(_totalSize)
        };

    UniquePtr<float[]> _patchScratch{makeUnique<float[]>(_totalSize)};
    UniquePtr<Complex<float>[]> _patchScratchOut {
        makeUnique<Complex<float>[]>(_totalSize)
        };

    UniquePtr<Complex<float>[]> _tmp { makeUnique<Complex<float>[]>(_totalSize) };
    UniquePtr<float[]> _tmpOut{makeUnique<float[]>(_totalSize)};
};
}  // namespace mc::dsp
