#pragma once

#include "lt/dsp/fft/FFT.hpp"

struct FFTConvolver {
    FFTConvolver(std::size_t signalSize, std::size_t patchSize);

    auto fft(double const* signal, double const* patch, double* output) const -> void;

    static auto direct(double const* signal, std::size_t n, double const* patch, std::size_t l, double* output) noexcept -> void;

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
