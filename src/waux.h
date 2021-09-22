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

auto mean(double const* vec, int n) -> double;

auto var(double const* vec, int n) -> double;

auto median(double* x, int n) -> double;

auto minindex(double const* arr, int n) -> int;

auto autocovar(double const* vec, int n, double* acov, int m) -> void;

auto autocorr(double const* vec, int n, double* acorr, int m) -> void;

#endif /* AUXILIARY_WAUX_H_ */
