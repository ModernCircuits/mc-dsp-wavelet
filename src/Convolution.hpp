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