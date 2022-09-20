// SPDX-License-Identifier: BSL-1.0

#include "PFFFT.hpp"

#include <mc/core/algorithm.hpp>
#include <mc/core/cassert.hpp>
#include <mc/core/cmath.hpp>
#include <mc/core/utility.hpp>

namespace mc::dsp {

PFFFT_Complex_Float::PFFFT_Complex_Float(std::size_t size)
    : _setup{pffft_new_setup(size, PFFFT_COMPLEX)}
{}

PFFFT_Complex_Float::~PFFFT_Complex_Float()
{
    if (_setup != nullptr) { pffft_destroy_setup(_setup); }
}

PFFFT_Complex_Float::PFFFT_Complex_Float(PFFFT_Complex_Float&& other) noexcept
    : _setup{std::exchange(other._setup, nullptr)}
{}

auto PFFFT_Complex_Float::operator=(PFFFT_Complex_Float&& other) noexcept
    -> PFFFT_Complex_Float&
{
    _setup = std::exchange(other._setup, nullptr);
    return *this;
}

auto PFFFT_Complex_Float::fft(Complex<float> const* in, Complex<float>* out) -> void
{
    auto const* input = reinterpret_cast<float const*>(in);  // NOLINT
    auto* output      = reinterpret_cast<float*>(out);       // NOLINT
    pffft_transform_ordered(_setup, input, output, nullptr, PFFFT_FORWARD);
}

auto PFFFT_Complex_Float::ifft(Complex<float> const* in, Complex<float>* out) -> void
{
    auto const* input = reinterpret_cast<float const*>(in);  // NOLINT
    auto* output      = reinterpret_cast<float*>(out);       // NOLINT
    pffft_transform_ordered(_setup, input, output, nullptr, PFFFT_BACKWARD);
}

PFFFT_Real_Float::PFFFT_Real_Float(size_t n)
    : _n{static_cast<int>(n)}
    , _setup{pffft_new_setup(static_cast<int>(n), PFFFT_REAL)}
{
    _tmp.resize(n);
}

PFFFT_Real_Float::~PFFFT_Real_Float()
{
    if (_setup != nullptr) { pffft_destroy_setup(_setup); }
}

PFFFT_Real_Float::PFFFT_Real_Float(PFFFT_Real_Float&& other) noexcept
    : _n{std::exchange(other._n, 0)}
    , _setup{std::exchange(other._setup, nullptr)}
    , _tmp{std::exchange(other._tmp, {})}
{}

auto PFFFT_Real_Float::operator=(PFFFT_Real_Float&& other) noexcept -> PFFFT_Real_Float&
{
    _n     = std::exchange(other._n, 0);
    _setup = std::exchange(other._setup, nullptr);
    _tmp   = std::exchange(other._tmp, {});
    return *this;
}

auto PFFFT_Real_Float::rfft(float const* inp, Complex<float>* oup) -> void
{
    pffft_transform_ordered(_setup, inp, (float*)oup, nullptr, PFFFT_FORWARD);

    // Move compressed DC/Nyquist components to correct location
    auto const h = _n / 2;
    oup[h]       = {oup[0].imag(), 0.0};
    oup[0]       = {oup[0].real(), 0.0};

    // Fill upper half with conjugate
    for (auto i = h + 1; i < _n; ++i) { oup[i] = std::conj(oup[_n - i]); }
}

auto PFFFT_Real_Float::irfft(Complex<float> const* inp, float* oup) -> void
{
    // Move DC/Nyquist components to compressed location
    std::copy(inp, inp + _n, std::begin(_tmp));
    _tmp[0] = {_tmp[0].real(), _tmp[_n / 2].real()};

    auto const* in = reinterpret_cast<float const*>(_tmp.data());
    pffft_transform_ordered(_setup, in, oup, nullptr, PFFFT_BACKWARD);
}

}  // namespace mc::dsp
