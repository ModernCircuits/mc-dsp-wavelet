#include "FFT.hpp"

#include <mc/core/cassert.hpp>
#include <mc/core/cmath.hpp>
#include <mc/core/format.hpp>
#include <mc/core/memory.hpp>

namespace {
constexpr auto pi2 = static_cast<float>(6.28318530717958647692528676655900577);
}  // namespace

namespace mc::dsp {

RFFT::RFFT(int n, FFTDirection direction)
{
    data_ = makeUnique<Complex<float>[]>(n / 2);
    cobj_ = makeUnique<FFT<float, KissFFT>>(n / 2, direction);

    for (auto k = 0; k < n / 2; ++k) {
        auto const theta = pi2 * k / n;
        data_[k].real(cos(theta));
        data_[k].imag(sin(theta));
    }
}

auto RFFT::performRealToComplex(float const* inp, Complex<float>* oup) -> void
{
    auto const n2 = cobj_->size();
    auto cinp     = makeUnique<Complex<float>[]>(n2);
    auto coup     = makeUnique<Complex<float>[]>(n2);

    for (auto i = 0; i < n2; ++i) {
        cinp[i].real(inp[2 * i]);
        cinp[i].imag(inp[2 * i + 1]);
    }

    cobj_->perform(cinp.get(), coup.get());

    oup[0].real(coup[0].real() + coup[0].imag());
    oup[0].imag(0.0F);

    for (auto i = 1; i < n2; ++i) {
        auto const temp1 = coup[i].imag() + coup[n2 - i].imag();
        auto const temp2 = coup[n2 - i].real() - coup[i].real();
        oup[i].real(
            (coup[i].real() + coup[n2 - i].real() + (temp1 * data_[i].real())
             + (temp2 * data_[i].imag()))
            / 2.0F
        );
        oup[i].imag(
            (coup[i].imag() - coup[n2 - i].imag() + (temp2 * data_[i].real())
             - (temp1 * data_[i].imag()))
            / 2.0F
        );
    }

    auto const n = n2 * 2;
    oup[n2].real(coup[0].real() - coup[0].imag());
    oup[n2].imag(0.0F);

    for (auto i = 1; i < n2; ++i) {
        oup[n - i].real(oup[i].real());
        oup[n - i].imag(-oup[i].imag());
    }
}

auto RFFT::performComplexToReal(Complex<float> const* inp, float* oup) -> void
{
    auto const n = static_cast<std::size_t>(cobj_->size());
    auto cinp    = makeUnique<Complex<float>[]>(n);
    auto coup    = makeUnique<Complex<float>[]>(n);

    for (auto i = std::size_t{0}; i < n; ++i) {
        auto const temp1 = -inp[i].imag() - inp[n - i].imag();
        auto const temp2 = -inp[n - i].real() + inp[i].real();
        cinp[i].real(
            inp[i].real() + inp[n - i].real() + (temp1 * data_[i].real())
            - (temp2 * data_[i].imag())
        );
        cinp[i].imag(
            inp[i].imag() - inp[n - i].imag() + (temp2 * data_[i].real())
            + (temp1 * data_[i].imag())
        );
    }

    cobj_->perform(cinp.get(), coup.get());

    for (auto i = std::size_t{0}; i < n; ++i) {
        oup[2 * i]     = coup[i].real();
        oup[2 * i + 1] = coup[i].imag();
    }
}
}  // namespace mc::dsp
