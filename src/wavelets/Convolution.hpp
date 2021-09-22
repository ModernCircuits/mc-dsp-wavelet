#pragma once

#include "FFT.hpp"
#include "real.h"

struct Convolution {
    std::unique_ptr<FftRealSet> fobj;
    std::unique_ptr<FftRealSet> iobj;
    int ilen1;
    int ilen2;
    int clen;
};

auto convInit(int n, int l) -> std::unique_ptr<Convolution>;

auto convDirect(fft_type const* inp1, int n, fft_type const* inp2, int l, fft_type* oup) -> void;
auto convFft(Convolution const& obj, fft_type const* inp1, fft_type const* inp2, fft_type* oup) -> void;