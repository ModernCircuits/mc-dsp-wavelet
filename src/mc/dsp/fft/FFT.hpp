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

template<typename T, typename BackendTag>
struct FFT
{
    using value_type   = T;
    using backend_type = typename FFTBackend<T, BackendTag>::handle_type;

    FFT(int n, FFTDirection direction);
    ~FFT();

    [[nodiscard]] auto direction() const noexcept -> FFTDirection;
    [[nodiscard]] auto size() const noexcept -> int;

    auto perform(Complex<T> const* inp, Complex<T>* output) -> void;

private:
    int _size;
    FFTDirection _direction;
    backend_type _engine;
};

template<typename T, typename BackendTag>
FFT<T, BackendTag>::FFT(int n, FFTDirection direction)
    : _size{n}
    , _direction{direction}
    , _engine{FFTBackend<T, BackendTag>::construct(static_cast<std::size_t>(n), direction)}
{}

template<typename T, typename BackendTag>
FFT<T, BackendTag>::~FFT()
{
    FFTBackend<T, BackendTag>::destroy(_engine);
}

template<typename T, typename BackendTag>
auto FFT<T, BackendTag>::direction() const noexcept -> FFTDirection
{
    return _direction;
}

template<typename T, typename BackendTag>
auto FFT<T, BackendTag>::size() const noexcept -> int
{
    return _size;
}

template<typename T, typename BackendTag>
auto FFT<T, BackendTag>::perform(Complex<T> const* input, Complex<T>* output) -> void
{
    FFTBackend<T, BackendTag>::perform(_engine, input, output, _direction);
}

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
