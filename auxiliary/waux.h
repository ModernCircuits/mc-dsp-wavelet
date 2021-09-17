/*
 * waux.h
 *
 *  Created on: Aug 22, 2017
 *      Author: Rafat Hussain
 */

#ifndef AUXILIARY_WAUX_H_
#define AUXILIARY_WAUX_H_

#include "wavelib.h"
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

auto compare_double(const void* a, const void* b) -> int;

auto mean(const double* vec, int N) -> double;

auto var(const double* vec, int N) -> double;

auto median(double* x, int N) -> double;

auto minindex(const double* arr, int N) -> int;

void getDWTAppx(wt_object wt, double* appx, int N);

void getDWTDetail(wt_object wt, double* detail, int N, int level);

void autocovar(const double* vec, int N, double* acov, int M);

void autocorr(const double* vec, int N, double* acorr, int M);

#endif /* AUXILIARY_WAUX_H_ */
