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

    auto rfft(float const* inp, Complex<float>* oup) -> void;
    auto irfft(Complex<float> const* inp, float* oup) -> void;

private:
    UniquePtr<FFT<float, KissFFT>> _cobj;
    UniquePtr<Complex<float>[]> _data;
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
