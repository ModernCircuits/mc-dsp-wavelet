/*
 * conv.h
 *
 *  Created on: May 1, 2013
 *      Author: Rafat Hussain
 */

#ifndef CONV_H_
#define CONV_H_

#include "real.h"

auto convDirect(fft_type const* inp1, int n, fft_type const* inp2, int l, fft_type* oup) -> void;
auto convFft(Convolution const& obj, fft_type const* inp1, fft_type const* inp2, fft_type* oup) -> void;

#endif /* CONV_H_ */
