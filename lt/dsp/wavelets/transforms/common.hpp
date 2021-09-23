#pragma once

#include "lt/dsp/wavelets/filters/utility.hpp"

auto dwtPerStride(double const* inp, int n, double const* lpd, double const* hpd, int lpdLen,
    double* cA, int lenCA, double* cD, int istride, int ostride) -> void;

auto dwtSymStride(double const* inp, int n, double const* lpd, double const* hpd, int lpdLen,
    double* cA, int lenCA, double* cD, int istride, int ostride) -> void;

auto modwtPerStride(int m, double const* inp, int n, double const* filt, int lpdLen,
    double* cA, int lenCA, double* cD, int istride, int ostride) -> void;

auto swtPerStride(int m, double const* inp, int n, double const* lpd, double const* hpd, int lpdLen,
    double* cA, int lenCA, double* cD, int istride, int ostride) -> void;

auto idwtPerStride(double const* cA, int lenCA, double const* cD, double const* lpr, double const* hpr,
    int lprLen, double* x, int istride, int ostride) -> void;

auto idwtSymStride(double const* cA, int lenCA, double const* cD, double const* lpr, double const* hpr,
    int lprLen, double* x, int istride, int ostride) -> void;

auto testSWTlength(int n, int j) -> int;

auto wmaxiter(int sigLen, int filtLen) -> int;
