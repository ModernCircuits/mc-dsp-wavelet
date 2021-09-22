#pragma once

#include <algorithm>
#include <memory>

#ifndef fft_type
#define fft_type double
#endif

#ifndef cplx_type
#define cplx_type double
#endif

template <typename T>
auto makeZeros(std::size_t length) -> std::unique_ptr<T[]>
{
    auto ptr = std::make_unique<T[]>(length);
    std::fill(ptr.get(), ptr.get() + length, T {});
    return ptr;
}

struct CplxData {
    cplx_type re;
    cplx_type im;
};

struct FftData {
    fft_type re;
    fft_type im;
};

struct FftSet {
    int N;
    int sgn;
    int factors[64];
    int lf;
    int lt;
    std::unique_ptr<FftData[]> data;
};

auto fftInit(int n, int sgn) -> std::unique_ptr<FftSet>;

struct FftRealSet {
    std::unique_ptr<FftSet> cobj;
    std::unique_ptr<FftData[]> data;
};

auto fftRealInit(int n, int sgn) -> std::unique_ptr<FftRealSet>;