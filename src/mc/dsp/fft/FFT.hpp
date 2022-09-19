// SPDX-License-Identifier: BSL-1.0

#pragma once

#include <mc/core/config.hpp>

#include <mc/dsp/fft/backend/KissFFT.hpp>
#include <mc/dsp/fft/backend/PFFFT.hpp>
#include <mc/dsp/fft/FFTBackend.hpp>
#include <mc/dsp/fft/FFTDirection.hpp>

#include <mc/core/algorithm.hpp>
#include <mc/core/complex.hpp>
#include <mc/core/memory.hpp>

namespace mc::dsp {

template<typename T, typename Engine>
struct FFT
{
    using value_type   = T;
    using backend_type = typename FFTBackend<T, Engine>::handle_type;

    FFT(int n);
    ~FFT();

    [[nodiscard]] auto size() const noexcept -> int;

    auto fft(Complex<T> const* inp, Complex<T>* output) -> void;
    auto ifft(Complex<T> const* inp, Complex<T>* output) -> void;

private:
    std::size_t _size;
    backend_type _fft;
    backend_type _ifft;
};

template<typename Engine, typename T>
inline auto fft(Engine& engine, Complex<T> const* input, Complex<T>* output)
    -> decltype(engine.fft(input, output))
{
    return engine.fft(input, output);
}

template<typename Engine, typename T>
inline auto ifft(Engine& engine, Complex<T> const* input, Complex<T>* output)
    -> decltype(engine.ifft(input, output))
{
    return engine.ifft(input, output);
}

template<typename T, typename Engine>
FFT<T, Engine>::FFT(int n)
    : _size{static_cast<std::size_t>(n)}
    , _fft{FFTBackend<T, Engine>::construct(_size, FFTDirection::forward)}
    , _ifft{FFTBackend<T, Engine>::construct(_size, FFTDirection::backward)}
{}

template<typename T, typename Engine>
FFT<T, Engine>::~FFT()
{
    FFTBackend<T, Engine>::destroy(_fft);
    FFTBackend<T, Engine>::destroy(_ifft);
}

template<typename T, typename Engine>
auto FFT<T, Engine>::size() const noexcept -> int
{
    return static_cast<int>(_size);
}

template<typename T, typename Engine>
auto FFT<T, Engine>::fft(Complex<T> const* input, Complex<T>* output) -> void
{
    FFTBackend<T, Engine>::perform(_fft, input, output, FFTDirection::forward);
}

template<typename T, typename Engine>
auto FFT<T, Engine>::ifft(Complex<T> const* input, Complex<T>* output) -> void
{
    FFTBackend<T, Engine>::perform(_ifft, input, output, FFTDirection::backward);
}

}  // namespace mc::dsp
