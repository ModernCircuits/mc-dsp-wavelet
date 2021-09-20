/*
 * conv.h
 *
 *  Created on: May 1, 2013
 *      Author: Rafat Hussain
 */

#ifndef CONV_H_
#define CONV_H_

#include "real.h"

auto factorf(int M) -> int;
auto findnext(int M) -> int;
auto findnexte(int M) -> int;

void conv_direct(fft_type const* inp1, int N, fft_type const* inp2, int L, fft_type* oup);
void conv_fft(conv_set const* obj, fft_type const* inp1, fft_type const* inp2, fft_type* oup);

#endif /* CONV_H_ */
