// SPDX-License-Identifier: BSL-1.0

#pragma once

#include <mc/core/config.hpp>

#include <mc/dsp/fft/FFT.hpp>

#include <mc/core/complex.hpp>
#include <mc/core/memory.hpp>
#include <mc/core/vector.hpp>

#include <pffft.h>

namespace mc::dsp {

struct RFFT
{
    RFFT(int n);
    ~RFFT();

    RFFT(RFFT const& other)                    = delete;
    auto operator=(RFFT const& other) -> RFFT& = delete;

    auto rfft(float const* inp, Complex<float>* oup) -> void;
    auto irfft(Complex<float> const* inp, float* oup) -> void;

private:
    int _n;
    PFFFT_Setup* _fft;

    // Needed, because we need to rearrange in inverse transform input
    // and it is passed as const.
    Vector<Complex<float>> _tmp;
};

template<typename Engine, typename T>
inline auto rfft(Engine& engine, T const* input, Complex<T>* output)
    -> decltype(engine.rfft(input, output))
{
    return engine.rfft(input, output);
}

template<typename Engine, typename T>
inline auto irfft(Engine& engine, Complex<T> const* input, T* output)
    -> decltype(engine.irfft(input, output))
{
    return engine.irfft(input, output);
}

}  // namespace mc::dsp
