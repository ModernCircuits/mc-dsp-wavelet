#include "FFTConvolver.hpp"

#include "lt/cmath.hpp"
#include <algorithm>
#include <cstdio>
#include <cstring>
#include <memory>

namespace {
[[nodiscard]] auto factorf(std::size_t m) -> std::size_t
{
    auto n = m;
    while (n % 7 == 0) {
        n = n / 7;
    }
    while (n % 3 == 0) {
        n = n / 3;
    }
    while (n % 5 == 0) {
        n = n / 5;
    }
    while (n % 2 == 0) {
        n = n / 2;
    }
    return n;
}

[[nodiscard]] auto findnexte(std::size_t m) -> std::size_t
{
    auto n = m;
    while (factorf(n) != 1 || n % 2 != 0) {
        ++n;
    }
    return n;
}
}

FFTConvolver::FFTConvolver(std::size_t signalSize, std::size_t patchSize)
    : signalSize_ { signalSize }
    , patchSize_ { patchSize }
    , totalSize_ { findnexte(signalSize + patchSize_ - 1U) }
    , forwardFFT_ { std::make_unique<RealFFT>(totalSize_, FFT::forward) }
    , backwardFFT_ { std::make_unique<RealFFT>(totalSize_, FFT::backward) }
{
}

auto FFTConvolver::convolute(double const* signal, double const* patch, double* output) const -> void
{
    std::fill(signalScratch_.get(), std::next(signalScratch_.get(), totalSize_), 0.0);
    std::fill(patchScratch_.get(), std::next(patchScratch_.get(), totalSize_), 0.0);
    std::fill(tmp_.get(), std::next(tmp_.get(), totalSize_), 0.0);
    std::fill(signalScratchOut_.get(), std::next(signalScratchOut_.get(), totalSize_), Complex<double> {});
    std::fill(patchScratchOut_.get(), std::next(patchScratchOut_.get(), totalSize_), Complex<double> {});
    std::fill(tmpOut_.get(), std::next(tmpOut_.get(), totalSize_), 0.0);

    for (auto i = std::size_t { 0 }; i < totalSize_; i++) {
        if (i < signalSize_) {
            signalScratch_[i] = signal[i];
        }
        if (i < patchSize_) {
            patchScratch_[i] = patch[i];
        }
    }

    forwardFFT_->performRealToComplex(signalScratch_.get(), signalScratchOut_.get());
    forwardFFT_->performRealToComplex(patchScratch_.get(), patchScratchOut_.get());

    for (auto i = std::size_t { 0 }; i < totalSize_; i++) {
        tmp_[i] = signalScratchOut_[i] * patchScratchOut_[i];
    }

    backwardFFT_->performComplexToReal(tmp_.get(), tmpOut_.get());

    auto const ls = signalSize_ + patchSize_ - 1U;
    for (auto i = std::size_t { 0 }; i < ls; i++) {
        output[i] = tmpOut_[i] / totalSize_;
    }
}
