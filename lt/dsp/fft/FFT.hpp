#pragma once

#include <algorithm>
#include <complex>
#include <memory>

#define PI2 6.28318530717958647692528676655900577

template <typename T>
auto makeZeros(std::size_t length) -> std::unique_ptr<T[]>
{
    auto ptr = std::make_unique<T[]>(length);
    std::fill(ptr.get(), ptr.get() + length, T {});
    return ptr;
}

template <typename T>
using Complex = std::complex<T>;

struct FFT {
    FFT(int n, int sgn);

    auto perform(Complex<double> const* inp, Complex<double>* oup) -> void;

    int N;
    int sgn;
    int factors[64];
    int lf;
    int lt;
    std::unique_ptr<Complex<double>[]> data;
};

struct FftRealSet {
    FftRealSet(int n, int sgn);

    auto performRealToComplex(double const* inp, Complex<double>* oup) -> void;
    auto performComplexToReal(Complex<double> const* inp, double* oup) -> void;

private:
    std::unique_ptr<FFT> cobj;
    std::unique_ptr<Complex<double>[]> data;
};

auto divideby(int m, int d) -> int;
