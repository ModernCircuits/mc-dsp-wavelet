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

#ifdef __cplusplus
extern "C" {
#endif

#define PI2 6.28318530717958647692528676655900577

#ifndef fft_type
#define fft_type double
#endif

/*
#define SADD(a,b) ((a)+(b))

#define SSUB(a,b) ((a)+(b))

#define SMUL(a,b) ((a)*(b))
*/

auto fft_init(int N, int sgn) -> fft_object;

void fft_exec(fft_object obj, fft_data* inp, fft_data* oup);

auto divideby(int M, int d) -> int;

auto dividebyN(int N) -> int;

//void arrrev(int M, int* arr);

auto factors(int M, int* arr) -> int;

void twiddle(fft_data* vec, int N, int radix);

void longvectorN(fft_data* sig, const int* array, int tx);

void free_fft(fft_object object);

#ifdef __cplusplus
}
#endif

#endif /* HSFFT_H_ */
