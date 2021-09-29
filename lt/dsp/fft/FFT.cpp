#include "FFT.hpp"

#include "lt/cassert.hpp"
#include "lt/cmath.hpp"
#include "lt/format.hpp"

#include "lt/memory.hpp"

namespace {

auto toKissFFTDirection(FFT::Direction direction) noexcept -> bool
{
    return direction != FFT::Direction::forward;
}

}

FFT::FFT(int n, Direction direction)
    : size_ { n }
    , direction_ { direction }
    , fftEngine_ { static_cast<std::size_t>(n), toKissFFTDirection(direction) }
{
}

auto FFT::direction() const noexcept -> Direction { return direction_; }

auto FFT::size() const noexcept -> int { return size_; }

auto FFT::perform(Complex<double> const* input, Complex<double>* ouput) -> void
{
    fftEngine_.transform(input, ouput);
}

auto divideby(int m, int d) -> int
{
    while (m % d == 0) {
        m = m / d;
    }
    if (m == 1) {
        return 1;
    }
    return 0;
}

RealFFT::RealFFT(int n, FFT::Direction direction)
{
    data_ = std::make_unique<Complex<double>[]>(n / 2);
    cobj_ = std::make_unique<FFT>(n / 2, direction);

    for (auto k = 0; k < n / 2; ++k) {
        auto const theta = pi2 * k / n;
        data_[k].real(std::cos(theta));
        data_[k].imag(std::sin(theta));
    }
}

auto RealFFT::performRealToComplex(double const* inp, Complex<double>* oup) -> void
{
    auto const n2 = cobj_->size();
    auto cinp = std::make_unique<Complex<double>[]>(n2);
    auto coup = std::make_unique<Complex<double>[]>(n2);

    for (auto i = 0; i < n2; ++i) {
        cinp[i].real(inp[2 * i]);
        cinp[i].imag(inp[2 * i + 1]);
    }

    cobj_->perform(cinp.get(), coup.get());

    oup[0].real(coup[0].real() + coup[0].imag());
    oup[0].imag(0.0);

    for (auto i = 1; i < n2; ++i) {
        auto const temp1 = coup[i].imag() + coup[n2 - i].imag();
        auto const temp2 = coup[n2 - i].real() - coup[i].real();
        oup[i].real((coup[i].real() + coup[n2 - i].real() + (temp1 * data_[i].real()) + (temp2 * data_[i].imag())) / 2.0);
        oup[i].imag((coup[i].imag() - coup[n2 - i].imag() + (temp2 * data_[i].real()) - (temp1 * data_[i].imag())) / 2.0);
    }

    auto const n = n2 * 2;
    oup[n2].real(coup[0].real() - coup[0].imag());
    oup[n2].imag(0.0);

    for (auto i = 1; i < n2; ++i) {
        oup[n - i].real(oup[i].real());
        oup[n - i].imag(-oup[i].imag());
    }
}

auto RealFFT::performComplexToReal(Complex<double> const* inp, double* oup) -> void
{
    auto const n = static_cast<std::size_t>(cobj_->size());
    auto cinp = std::make_unique<Complex<double>[]>(n);
    auto coup = std::make_unique<Complex<double>[]>(n);

    for (auto i = std::size_t { 0 }; i < n; ++i) {
        auto const temp1 = -inp[i].imag() - inp[n - i].imag();
        auto const temp2 = -inp[n - i].real() + inp[i].real();
        cinp[i].real(inp[i].real() + inp[n - i].real() + (temp1 * data_[i].real()) - (temp2 * data_[i].imag()));
        cinp[i].imag(inp[i].imag() - inp[n - i].imag() + (temp2 * data_[i].real()) + (temp1 * data_[i].imag()));
    }

    cobj_->perform(cinp.get(), coup.get());

    for (auto i = std::size_t { 0 }; i < n; ++i) {
        oup[2 * i] = coup[i].real();
        oup[2 * i + 1] = coup[i].imag();
    }
}
