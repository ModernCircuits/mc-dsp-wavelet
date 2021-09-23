#pragma once

#include "lt/dsp/convolution/convolute.hpp"
#include "lt/dsp/fft/FFT.hpp"

struct FFTConvolver {
    using value_type = double;

    FFTConvolver(std::size_t signalSize, std::size_t patchSize);

    auto convolute(double const* signal, double const* patch, double* output) const -> void;

private:
    std::size_t signalSize_;
    std::size_t patchSize_;
    std::size_t totalSize_;

    std::unique_ptr<RealFFT> forwardFFT_;
    std::unique_ptr<RealFFT> backwardFFT_;

    std::unique_ptr<double[]> signalScratch_ { std::make_unique<double[]>(totalSize_) };
    std::unique_ptr<Complex<double>[]> signalScratchOut_ { std::make_unique<Complex<double>[]>(totalSize_) };

    std::unique_ptr<double[]> patchScratch_ { std::make_unique<double[]>(totalSize_) };
    std::unique_ptr<Complex<double>[]> patchScratchOut_ { std::make_unique<Complex<double>[]>(totalSize_) };

    std::unique_ptr<Complex<double>[]> tmp_ { std::make_unique<Complex<double>[]>(totalSize_) };
    std::unique_ptr<double[]> tmpOut_ { std::make_unique<double[]>(totalSize_) };
};
