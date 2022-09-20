// SPDX-License-Identifier: BSL-1.0

#include "RFFT.hpp"

#include <mc/core/cassert.hpp>
#include <mc/core/cmath.hpp>
#include <mc/core/format.hpp>

namespace mc::dsp {

RFFT::RFFT(int n) : _n{n}, _fft{pffft_new_setup(n, PFFFT_REAL)} { _tmp.resize(n); }

RFFT::~RFFT()
{
    if (_fft != nullptr) { pffft_destroy_setup(_fft); }
}

auto RFFT::rfft(float const* inp, Complex<float>* oup) -> void
{
    pffft_transform_ordered(_fft, inp, (float*)oup, nullptr, PFFFT_FORWARD);

    // Move compressed DC/Nyquist components to correct location
    auto const h = _n / 2;
    oup[h]       = {oup[0].imag(), 0.0};
    oup[0]       = {oup[0].real(), 0.0};

    // Fill upper half with conjugate
    for (auto i = h + 1; i < _n; ++i) { oup[i] = std::conj(oup[_n - i]); }
}

auto RFFT::irfft(Complex<float> const* inp, float* oup) -> void
{
    // Move DC/Nyquist components to compressed location
    ranges::copy(inp, inp + _n, ranges::begin(_tmp));
    _tmp[0] = {_tmp[0].real(), _tmp[_n / 2].real()};

    pffft_transform_ordered(_fft, (float const*)_tmp.data(), oup, nullptr, PFFFT_BACKWARD);
}
}  // namespace mc::dsp
