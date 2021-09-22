/*
 * hsfft.h
 *
 *  Created on: Apr 14, 2013
 *      Author: Rafat Hussain
 */

#ifndef HSFFT_H_
#define HSFFT_H_

#include <cmath>
#include <cstdlib>
#include <cstring>

#include "wavelib.h"

#define PI2 6.28318530717958647692528676655900577

#ifndef fft_type
#define fft_type double
#endif

auto fftExec(FftSet& obj, FftData* inp, FftData* oup) -> void;

auto divideby(int m, int d) -> int;

auto dividebyN(int n) -> int;

auto factors(int m, int* arr) -> int;

auto twiddle(FftData* vec, int n, int radix) -> void;

auto longvectorN(FftData* sig, int const* array, int tx) -> void;

#endif /* HSFFT_H_ */
