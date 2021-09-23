#pragma once

#include <algorithm>
#include <memory>

#define PI2 6.28318530717958647692528676655900577

template <typename T>
auto makeZeros(std::size_t length) -> std::unique_ptr<T[]>
{
    auto ptr = std::make_unique<T[]>(length);
    std::fill(ptr.get(), ptr.get() + length, T {});
    return ptr;
}

struct Complex {
    double re;
    double im;
};

struct FFT {
    FFT(int n, int sgn);

    auto perform(Complex* inp, Complex* oup) -> void;

    int N;
    int sgn;
    int factors[64];
    int lf;
    int lt;
    std::unique_ptr<Complex[]> data;
};

struct FftRealSet {
    std::unique_ptr<FFT> cobj;
    std::unique_ptr<Complex[]> data;
};

auto fftRealInit(int n, int sgn) -> std::unique_ptr<FftRealSet>;

auto fftR2cExec(FftRealSet* obj, double const* inp, Complex* oup) -> void;
auto fftC2rExec(FftRealSet* obj, Complex* inp, double* oup) -> void;

auto divideby(int m, int d) -> int;
