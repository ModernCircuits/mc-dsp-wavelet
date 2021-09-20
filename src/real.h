/*
 * real.h
 *
 *  Created on: Apr 20, 2013
 *      Author: Rafat Hussain
 */

#ifndef REAL_H_
#define REAL_H_

#include "hsfft.h"

void fft_r2c_exec(fft_real_set* obj, fft_type const* inp, fft_data* oup);
void fft_c2r_exec(fft_real_set* obj, fft_data* inp, fft_type* oup);
void free_real_fft(fft_real_set* object);

#endif /* REAL_H_ */
