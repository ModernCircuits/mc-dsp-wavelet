#pragma once

#include "lt/dsp/fft/FFT.hpp"

struct Convolution {
    Convolution(int n, int l);

    std::unique_ptr<RealFFT> fobj;
    std::unique_ptr<RealFFT> iobj;
    int ilen1;
    int ilen2;
    int clen;
};

auto convDirect(double const* inp1, int n, double const* inp2, int l, double* oup) -> void;
auto convFft(Convolution const& obj, double const* inp1, double const* inp2, double* oup) -> void;
