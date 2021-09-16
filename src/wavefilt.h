/*
Copyright (c) 2014, Rafat Hussain
Copyright (c) 2016, Holger Nahrstaedt
*/
#ifndef WAVEFILT_H_
#define WAVEFILT_H_

#include "conv.h"
#include <cstdio>
#define USE_MATH_DEFINES
#include <cmath>

#ifdef __cplusplus
extern "C" {
#endif

auto filtlength(const char* name) -> int;

auto filtcoef(const char* name, double* lp1, double* hp1, double* lp2, double* hp2) -> int;

void copy_reverse(const double* in, int N, double* out);
void qmf_even(const double* in, int N, double* out);
void qmf_wrev(const double* in, int N, double* out);
void copy(const double* in, int N, double* out);
#ifdef __cplusplus
}
#endif

#endif /* WAVEFILT_H_ */