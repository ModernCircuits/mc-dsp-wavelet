// SPDX-License-Identifier: BSL-1.0

#pragma once

#include <mc/core/config.hpp>

#include <mc/dsp/fft/FFT.hpp>

#include <mc/core/complex.hpp>
#include <mc/core/memory.hpp>

namespace mc::dsp {

struct RFFT
{
    RFFT(int n, FFTDirection direction);

    auto performRealToComplex(float const* inp, Complex<float>* oup) -> void;
    auto performComplexToReal(Complex<float> const* inp, float* oup) -> void;

private:
    UniquePtr<FFT<float, KissFFT>> _cobj;
    UniquePtr<Complex<float>[]> _data;
};

}  // namespace mc::dsp
