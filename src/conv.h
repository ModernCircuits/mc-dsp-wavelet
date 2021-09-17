/*
 * conv.h
 *
 *  Created on: May 1, 2013
 *      Author: Rafat Hussain
 */

#ifndef CONV_H_
#define CONV_H_

#include "real.h"

#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))

auto conv_init(int N, int L) -> conv_object;

auto factorf(int M) -> int;

auto findnext(int M) -> int;

auto findnexte(int M) -> int;

void conv_direct(const fft_type* inp1, int N, const fft_type* inp2, int L, fft_type* oup);

void conv_directx(const fft_type* inp1, int N, const fft_type* inp2, int L, fft_type* oup);

//void conv_fft(const conv_object obj,fft_type *inp1,fft_type *inp2,fft_type *oup);

//void conv_fft(const conv_object obj,fft_type *inp1,fft_type *inp2,fft_type *oup);

void conv_fft(conv_object obj, const fft_type* inp1, const fft_type* inp2, fft_type* oup);

//void free_conv(conv_object object);

void free_conv(conv_object object);

#endif /* CONV_H_ */
