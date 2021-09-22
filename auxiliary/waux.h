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

auto mean(double const* vec, int N) -> double;

auto var(double const* vec, int N) -> double;

auto median(double* x, int N) -> double;

auto minindex(double const* arr, int N) -> int;

auto autocovar(double const* vec, int N, double* acov, int M) -> void;

auto autocorr(double const* vec, int N, double* acorr, int M) -> void;

#endif /* AUXILIARY_WAUX_H_ */
