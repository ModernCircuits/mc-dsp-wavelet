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

void fft_exec(fft_object obj, fft_data* inp, fft_data* oup);

auto divideby(int M, int d) -> int;

auto dividebyN(int N) -> int;

auto factors(int M, int* arr) -> int;

void twiddle(fft_data* vec, int N, int radix);

void longvectorN(fft_data* sig, int const* array, int tx);

void free_fft(fft_object object);

#endif /* HSFFT_H_ */
