/*
 * real.h
 *
 *  Created on: Apr 20, 2013
 *      Author: Rafat Hussain
 */

#ifndef REAL_H_
#define REAL_H_

#include "hsfft.h"

auto fftR2cExec(FftRealSet* obj, fft_type const* inp, FftData* oup) -> void;
auto fftC2rExec(FftRealSet* obj, FftData* inp, fft_type* oup) -> void;

#endif /* REAL_H_ */
