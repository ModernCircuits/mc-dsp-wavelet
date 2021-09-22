/// \copyright Copyright (c) 2014, Rafat Hussain
/// \copyright Copyright (c) 2016, Holger Nahrstaedt

#ifndef WAVEFILT_H_
#define WAVEFILT_H_

#include "conv.h"
#include <cstdio>
#define USE_MATH_DEFINES
#include <cmath>

auto filtlength(char const* name) -> int;
auto filtcoef(char const* name, double* lp1, double* hp1, double* lp2, double* hp2) -> int;

void copyReverse(double const* in, int n, double* out);
void qmfEven(double const* in, int n, double* out);
void qmfWrev(double const* in, int n, double* out);
void copy(double const* in, int n, double* out);

#endif /* WAVEFILT_H_ */