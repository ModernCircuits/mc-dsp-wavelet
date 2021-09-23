#include "Convolution.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <memory>

namespace {
[[nodiscard]] auto factorf(int m) -> int
{
    auto n = m;
    while (n % 7 == 0) {
        n = n / 7;
    }
    while (n % 3 == 0) {
        n = n / 3;
    }
    while (n % 5 == 0) {
        n = n / 5;
    }
    while (n % 2 == 0) {
        n = n / 2;
    }
    return n;
}

[[nodiscard]] auto findnexte(int m) -> int
{
    auto n = m;
    while (factorf(n) != 1 || n % 2 != 0) {
        ++n;
    }
    return n;
}
}

Convolution::Convolution(int n, int l)
    : ilen1_ { n }
    , ilen2_ { l }
    , clen_ { findnexte(n + l - 1) }
    , fobj_ { std::make_unique<RealFFT>(clen_, 1) }
    , iobj_ { std::make_unique<RealFFT>(clen_, -1) }
{
}

auto Convolution::fft(double const* inp1, double const* inp2, double* oup) const -> void
{

    auto n = clen_;
    auto l1 = ilen1_;
    auto l2 = ilen2_;
    auto ls = l1 + l2 - 1;

    // auto a = std::make_unique<double[]>(n);
    // auto b = std::make_unique<double[]>(n);
    // auto c = std::make_unique<Complex<double>[]>(n);
    // auto ao = std::make_unique<Complex<double>[]>(n);
    // auto bo = std::make_unique<Complex<double>[]>(n);
    // auto co = std::make_unique<double[]>(n);

    std::fill(a.get(), std::next(a.get(), clen_), 0.0);
    std::fill(b.get(), std::next(b.get(), clen_), 0.0);
    std::fill(c.get(), std::next(c.get(), clen_), 0.0);
    std::fill(ao.get(), std::next(ao.get(), clen_), Complex<double> {});
    std::fill(bo.get(), std::next(bo.get(), clen_), Complex<double> {});
    std::fill(co.get(), std::next(co.get(), clen_), 0.0);

    for (auto i = 0; i < n; i++) {
        if (i < l1) {
            a[i] = inp1[i];
        } else {
            a[i] = 0.0;
        }

        if (i < l2) {
            b[i] = inp2[i];
        } else {
            b[i] = 0.0;
        }
    }

    fobj_->performRealToComplex(a.get(), ao.get());
    fobj_->performRealToComplex(b.get(), bo.get());

    for (auto i = 0; i < n; i++) {
        c[i] = Complex<double> {
            ao[i].real() * bo[i].real() - ao[i].imag() * bo[i].imag(),
            ao[i].imag() * bo[i].real() + ao[i].real() * bo[i].imag(),
        };
    }

    iobj_->performComplexToReal(c.get(), co.get());

    for (auto i = 0; i < ls; i++) {
        oup[i] = co[i] / n;
    }
}

auto Convolution::direct(double const* inp1, int n, double const* inp2, int l, double* oup) noexcept -> void
{
    auto const mm = n + l - 1;

    if (n >= l) {
        auto i = 0;
        for (auto k = 0; k < l; k++) {
            oup[k] = 0.0;
            for (auto m = 0; m <= k; m++) {
                oup[k] += inp1[m] * inp2[k - m];
            }
        }
        for (auto k = l; k < mm; k++) {
            oup[k] = 0.0;
            i++;
            auto const t1 = l + i;
            auto const tmin = std::min<double>(t1, n);
            for (auto m = i; m < tmin; m++) {
                oup[k] += inp1[m] * inp2[k - m];
            }
        }
        return;
    }

    auto i = 0;
    for (auto k = 0; k < n; k++) {
        oup[k] = 0.0;
        for (auto m = 0; m <= k; m++) {
            oup[k] += inp2[m] * inp1[k - m];
        }
    }
    for (auto k = n; k < mm; k++) {
        oup[k] = 0.0;
        i++;
        auto const t1 = n + i;
        auto const tmin = std::min<double>(t1, l);
        for (auto m = i; m < tmin; m++) {
            oup[k] += inp2[m] * inp1[k - m];
        }
    }
}