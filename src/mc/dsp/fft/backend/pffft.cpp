// SPDX-License-Identifier: BSL-1.0

#include "pffft.hpp"

#include <mc/core/algorithm.hpp>
#include <mc/core/cassert.hpp>
#include <mc/core/cmath.hpp>
#include <mc/core/utility.hpp>

namespace mc {

PFFFT_Complex_Float::PFFFT_Complex_Float(size_t size)
    : _setup{PFFFT_Handle{pffft_new_setup(static_cast<int>(size), PFFFT_COMPLEX)}}
{}

auto PFFFT_Complex_Float::fft(Span<Complex<float> const> in, Span<Complex<float>> out)
    -> void
{
    auto const* input = reinterpret_cast<float const*>(in.data());  // NOLINT
    auto* output      = reinterpret_cast<float*>(out.data());       // NOLINT
    pffft_transform_ordered(_setup.get(), input, output, nullptr, PFFFT_FORWARD);
}

auto PFFFT_Complex_Float::ifft(Span<Complex<float> const> in, Span<Complex<float>> out)
    -> void
{
    auto const* input = reinterpret_cast<float const*>(in.data());  // NOLINT
    auto* output      = reinterpret_cast<float*>(out.data());       // NOLINT
    pffft_transform_ordered(_setup.get(), input, output, nullptr, PFFFT_BACKWARD);
}

PFFFT_Real_Float::PFFFT_Real_Float(size_t n)
    : _n{static_cast<int>(n)}
    , _setup{PFFFT_Handle{pffft_new_setup(static_cast<int>(n), PFFFT_REAL)}}
{
    _tmp.resize(n);
}

auto PFFFT_Real_Float::rfft(Span<float const> inp, Span<Complex<float>> oup) -> void
{
    pffft_transform_ordered(
        _setup.get(),
        inp.data(),
        (float*)oup.data(),
        nullptr,
        PFFFT_FORWARD
    );

    // Move compressed DC/Nyquist components to correct location
    auto const h = _n / 2;
    oup[h]       = {oup[0].imag(), 0.0};
    oup[0]       = {oup[0].real(), 0.0};

    // Fill upper half with conjugate
    for (auto i = h + 1; i < _n; ++i) { oup[i] = std::conj(oup[_n - i]); }
}

auto PFFFT_Real_Float::irfft(Span<Complex<float> const> inp, Span<float> oup) -> void
{
    // Move DC/Nyquist components to compressed location
    ranges::copy(inp, ranges::begin(_tmp));
    _tmp[0] = {_tmp[0].real(), _tmp[_n / 2].real()};

    auto const* in = reinterpret_cast<float const*>(_tmp.data());
    pffft_transform_ordered(_setup.get(), in, oup.data(), nullptr, PFFFT_BACKWARD);
}

}  // namespace mc
