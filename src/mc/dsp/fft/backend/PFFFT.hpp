// SPDX-License-Identifier: BSL-1.0

#pragma once

#include <mc/dsp/fft/FFTDirection.hpp>

#include <mc/core/algorithm.hpp>
#include <mc/core/cassert.hpp>
#include <mc/core/cmath.hpp>
#include <mc/core/complex.hpp>
#include <mc/core/utility.hpp>
#include <mc/core/vector.hpp>

#include <pffft.h>

namespace mc::dsp {

// struct PFFFT
// {};

// template<>
// struct FFTBackend<float, PFFFT>
// {
//     using value_type  = float;
//     using handle_type = PFFFT_Setup*;

// static auto construct(std::size_t size, FFTDirection /*direction*/)
// {
//     return pffft_new_setup(static_cast<int>(size), PFFFT_COMPLEX);
// }

// static auto destroy(handle_type handle) { pffft_destroy_setup(handle); }

// static auto perform(
//     handle_type handle,
//     Complex<float> const* in,
//     Complex<float>* out,
//     FFTDirection direction
// ) -> void
// {
//     auto const dir
//         = direction == FFTDirection::forward ? PFFFT_FORWARD : PFFFT_BACKWARD;
//     auto const* input = reinterpret_cast<float const*>(in);  // NOLINT
//     auto* output      = reinterpret_cast<float*>(out);       // NOLINT
//     pffft_transform_ordered(handle, input, output, nullptr, dir);
// }
// };

struct PFFFT
{
    PFFFT(int size) : _setup{pffft_new_setup(size, PFFFT_COMPLEX)} {}

    ~PFFFT()
    {
        if (_setup != nullptr) { pffft_destroy_setup(_setup); }
    }

    PFFFT(PFFFT const& other)                    = delete;
    auto operator=(PFFFT const& other) -> PFFFT& = delete;

    PFFFT(PFFFT&& other) noexcept : _setup{std::exchange(other._setup, nullptr)} {}

    auto operator=(PFFFT&& other) noexcept -> PFFFT&
    {
        _setup = std::exchange(other._setup, nullptr);
        return *this;
    }

    auto fft(std::complex<float> const* in, std::complex<float>* out)
    {
        auto const* input = reinterpret_cast<float const*>(in);  // NOLINT
        auto* output      = reinterpret_cast<float*>(out);       // NOLINT
        pffft_transform_ordered(_setup, input, output, nullptr, PFFFT_FORWARD);
    }

    auto ifft(std::complex<float> const* in, std::complex<float>* out)
    {
        auto const* input = reinterpret_cast<float const*>(in);  // NOLINT
        auto* output      = reinterpret_cast<float*>(out);       // NOLINT
        pffft_transform_ordered(_setup, input, output, nullptr, PFFFT_BACKWARD);
    }

private:
    PFFFT_Setup* _setup;
};

struct PFFFT_Real
{
    PFFFT_Real(size_t n);
    ~PFFFT_Real();

    PFFFT_Real(PFFFT_Real const& other)                    = delete;
    auto operator=(PFFFT_Real const& other) -> PFFFT_Real& = delete;

    PFFFT_Real(PFFFT_Real&& other) noexcept;
    auto operator=(PFFFT_Real&& other) noexcept -> PFFFT_Real&;

    auto rfft(float const* inp, Complex<float>* oup) -> void;
    auto irfft(Complex<float> const* inp, float* oup) -> void;

private:
    int _n;
    PFFFT_Setup* _setup;

    // Needed, because we need to rearrange in inverse transform input
    // and it is passed as const.
    Vector<Complex<float>> _tmp;
};

inline PFFFT_Real::PFFFT_Real(size_t n)
    : _n{static_cast<int>(n)}
    , _setup{pffft_new_setup(static_cast<int>(n), PFFFT_REAL)}
{
    _tmp.resize(n);
}

inline PFFFT_Real::~PFFFT_Real()
{
    if (_setup != nullptr) { pffft_destroy_setup(_setup); }
}

inline PFFFT_Real::PFFFT_Real(PFFFT_Real&& other) noexcept
    : _n{std::exchange(other._n, 0)}
    , _setup{std::exchange(other._setup, nullptr)}
    , _tmp{std::exchange(other._tmp, {})}
{}

inline auto PFFFT_Real::operator=(PFFFT_Real&& other) noexcept -> PFFFT_Real&
{
    _n     = std::exchange(other._n, 0);
    _setup = std::exchange(other._setup, nullptr);
    _tmp   = std::exchange(other._tmp, {});
    return *this;
}

inline auto PFFFT_Real::rfft(float const* inp, Complex<float>* oup) -> void
{
    pffft_transform_ordered(_setup, inp, (float*)oup, nullptr, PFFFT_FORWARD);

    // Move compressed DC/Nyquist components to correct location
    auto const h = _n / 2;
    oup[h]       = {oup[0].imag(), 0.0};
    oup[0]       = {oup[0].real(), 0.0};

    // Fill upper half with conjugate
    for (auto i = h + 1; i < _n; ++i) { oup[i] = std::conj(oup[_n - i]); }
}

inline auto PFFFT_Real::irfft(Complex<float> const* inp, float* oup) -> void
{
    // Move DC/Nyquist components to compressed location
    std::copy(inp, inp + _n, std::begin(_tmp));
    _tmp[0] = {_tmp[0].real(), _tmp[_n / 2].real()};

    pffft_transform_ordered(
        _setup,
        (float const*)_tmp.data(),
        oup,
        nullptr,
        PFFFT_BACKWARD
    );
}

}  // namespace mc::dsp
