/*
 * real.h
 *
 *  Created on: Apr 20, 2013
 *      Author: Rafat Hussain
 */

#ifndef REAL_H_
#define REAL_H_

#include "hsfft.h"

auto fft_real_init(int N, int sgn) -> fft_real_object;
void fft_r2c_exec(fft_real_object obj, const fft_type* inp, fft_data* oup);
void fft_c2r_exec(fft_real_object obj, fft_data* inp, fft_type* oup);
void free_real_fft(fft_real_object object);

#endif /* REAL_H_ */
