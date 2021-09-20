/// \copyright Copyright (c) 2014, Rafat Hussain

#ifndef WTMATH_H_
#define WTMATH_H_

#include "wavefilt.h"

void dwt_per_stride(double const* inp, int N, double const* lpd, double const* hpd, int lpd_len,
    double* cA, int len_cA, double* cD, int istride, int ostride);

void dwt_sym_stride(double const* inp, int N, double const* lpd, double const* hpd, int lpd_len,
    double* cA, int len_cA, double* cD, int istride, int ostride);

void modwt_per_stride(int M, double const* inp, int N, double const* filt, int lpd_len,
    double* cA, int len_cA, double* cD, int istride, int ostride);

void swt_per_stride(int M, double const* inp, int N, double const* lpd, double const* hpd, int lpd_len,
    double* cA, int len_cA, double* cD, int istride, int ostride);

void idwt_per_stride(double const* cA, int len_cA, double const* cD, double const* lpr, double const* hpr,
    int lpr_len, double* X, int istride, int ostride);

void idwt_sym_stride(double const* cA, int len_cA, double const* cD, double const* lpr, double const* hpr,
    int lpr_len, double* X, int istride, int ostride);

void imodwt_per_stride(int M, double const* cA, int len_cA, double const* cD, double const* filt,
    int lf, double* X, int istride, int ostride);

void idwt2_shift(int shift, int rows, int cols, double* lpr, double* hpr, int lf,
    double* A, double* H, double* V, double* D, double* oup);

auto upsamp(double const* x, int lenx, int M, double* y) -> int;

auto upsamp2(double const* x, int lenx, int M, double* y) -> int;

auto downsamp(double const* x, int lenx, int M, double* y) -> int;

auto per_ext(double const* sig, int len, int a, double* oup) -> int;

auto symm_ext(double const* sig, int len, int a, double* oup) -> int;

void circshift(double* array, int N, int L);

auto testSWTlength(int N, int J) -> int;

auto wmaxiter(int sig_len, int filt_len) -> int;

auto costfunc(double* x, int N, char const* entropy, double p) -> double;

#endif /* WAVELIB_H_ */
