/*
 * real.h
 *
 *  Created on: Apr 20, 2013
 *      Author: Rafat Hussain
 */

#ifndef REAL_H_
#define REAL_H_

#include "hsfft.h"

void fftR2cExec(FftRealSet* obj, fft_type const* inp, FftData* oup);
void fftC2rExec(FftRealSet* obj, FftData* inp, fft_type* oup);

#endif /* REAL_H_ */
