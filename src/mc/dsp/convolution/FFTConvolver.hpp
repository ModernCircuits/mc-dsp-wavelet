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
    std::size_t signalSize_;
    std::size_t patchSize_;
    std::size_t totalSize_;

    UniquePtr<RFFT> forwardFFT_;
    UniquePtr<RFFT> backwardFFT_;

    UniquePtr<float[]> signalScratch_{makeUnique<float[]>(totalSize_)};
    UniquePtr<Complex<float>[]> signalScratchOut_ {
        makeUnique<Complex<float>[]>(totalSize_)
        };

    UniquePtr<float[]> patchScratch_{makeUnique<float[]>(totalSize_)};
    UniquePtr<Complex<float>[]> patchScratchOut_ {
        makeUnique<Complex<float>[]>(totalSize_)
        };

    UniquePtr<Complex<float>[]> tmp_ { makeUnique<Complex<float>[]>(totalSize_) };
    UniquePtr<float[]> tmpOut_{makeUnique<float[]>(totalSize_)};
};
}  // namespace mc::dsp
