// SPDX-License-Identifier: BSL-1.0

#pragma once

#include <mc/core/complex.hpp>
#include <mc/core/vector.hpp>

#include <pffft.h>

namespace mc::dsp {

struct PFFFT_Complex_Float
{
    explicit PFFFT_Complex_Float(std::size_t size);
    ~PFFFT_Complex_Float();

    PFFFT_Complex_Float(PFFFT_Complex_Float const& other)                    = delete;
    auto operator=(PFFFT_Complex_Float const& other) -> PFFFT_Complex_Float& = delete;

    PFFFT_Complex_Float(PFFFT_Complex_Float&& other) noexcept;
    auto operator=(PFFFT_Complex_Float&& other) noexcept -> PFFFT_Complex_Float&;

    auto fft(Complex<float> const* in, Complex<float>* out) -> void;
    auto ifft(Complex<float> const* in, Complex<float>* out) -> void;

private:
    PFFFT_Setup* _setup;
};

struct PFFFT_Real_Float
{
    explicit PFFFT_Real_Float(size_t n);
    ~PFFFT_Real_Float();

    PFFFT_Real_Float(PFFFT_Real_Float const& other)                    = delete;
    auto operator=(PFFFT_Real_Float const& other) -> PFFFT_Real_Float& = delete;

    PFFFT_Real_Float(PFFFT_Real_Float&& other) noexcept;
    auto operator=(PFFFT_Real_Float&& other) noexcept -> PFFFT_Real_Float&;

    auto rfft(float const* inp, Complex<float>* oup) -> void;
    auto irfft(Complex<float> const* inp, float* oup) -> void;

private:
    int _n;
    PFFFT_Setup* _setup;

    // We need to rearrange the inverse transform input
    // and it is passed as const. So we make a copy and modify.
    Vector<Complex<float>> _tmp;
};

}  // namespace mc::dsp
