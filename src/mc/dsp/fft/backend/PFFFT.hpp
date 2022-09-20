// SPDX-License-Identifier: BSL-1.0

#pragma once

#include <mc/core/complex.hpp>
#include <mc/core/memory.hpp>
#include <mc/core/vector.hpp>

#include <pffft.h>

namespace mc::dsp {

struct PFFFT_Deleter
{
    auto operator()(PFFFT_Setup* setup) { pffft_destroy_setup(setup); }
};

using PFFFT_Handle = UniquePtr<PFFFT_Setup, PFFFT_Deleter>;

struct PFFFT_Complex_Float
{
    explicit PFFFT_Complex_Float(std::size_t size);

    auto fft(Complex<float> const* in, Complex<float>* out) -> void;
    auto ifft(Complex<float> const* in, Complex<float>* out) -> void;

private:
    PFFFT_Handle _setup;
};

struct PFFFT_Real_Float
{
    explicit PFFFT_Real_Float(size_t n);

    auto rfft(float const* inp, Complex<float>* oup) -> void;
    auto irfft(Complex<float> const* inp, float* oup) -> void;

private:
    int _n;
    PFFFT_Handle _setup;

    // We need to rearrange the inverse transform input
    // and it is passed as const. So we make a copy and modify.
    Vector<Complex<float>> _tmp;
};

}  // namespace mc::dsp
