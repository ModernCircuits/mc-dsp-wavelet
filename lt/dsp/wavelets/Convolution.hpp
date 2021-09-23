#pragma once

#include "FFT.hpp"

struct Convolution {
    std::unique_ptr<FftRealSet> fobj;
    std::unique_ptr<FftRealSet> iobj;
    int ilen1;
    int ilen2;
    int clen;
};

auto convInit(int n, int l) -> std::unique_ptr<Convolution>;

auto convDirect(double const* inp1, int n, double const* inp2, int l, double* oup) -> void;
auto convFft(Convolution const& obj, double const* inp1, double const* inp2, double* oup) -> void;
