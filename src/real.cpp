/*
 * real.c
 *
 *  Created on: Apr 20, 2013
 *      Author: Rafat Hussain
 */
#include "real.h"

#include <cstdio>
#include <memory>

auto fftRealInit(int n, int sgn) -> std::unique_ptr<FftRealSet>
{
    auto obj = std::make_unique<FftRealSet>();
    obj->data = std::make_unique<FftData[]>(n / 2);
    obj->cobj = fftInit(n / 2, sgn);

    for (auto k = 0; k < n / 2; ++k) {
        auto const theta = PI2 * k / n;
        obj->data[k].re = cos(theta);
        obj->data[k].im = sin(theta);
    }

    return obj;
}

auto fftR2cExec(FftRealSet* obj, fft_type const* inp, FftData* oup) -> void
{
    int i;
    int n2;
    int n;
    fft_type temp1;
    fft_type temp2;
    n2 = obj->cobj->N;
    n = n2 * 2;

    auto cinp = std::make_unique<FftData[]>(n2);
    auto coup = std::make_unique<FftData[]>(n2);

    for (i = 0; i < n2; ++i) {
        cinp[i].re = inp[2 * i];
        cinp[i].im = inp[2 * i + 1];
    }

    fftExec(*(obj->cobj), cinp.get(), coup.get());

    oup[0].re = coup[0].re + coup[0].im;
    oup[0].im = 0.0;

    for (i = 1; i < n2; ++i) {
        temp1 = coup[i].im + coup[n2 - i].im;
        temp2 = coup[n2 - i].re - coup[i].re;
        oup[i].re = (coup[i].re + coup[n2 - i].re + (temp1 * obj->data[i].re) + (temp2 * obj->data[i].im)) / 2.0;
        oup[i].im = (coup[i].im - coup[n2 - i].im + (temp2 * obj->data[i].re) - (temp1 * obj->data[i].im)) / 2.0;
    }

    oup[n2].re = coup[0].re - coup[0].im;
    oup[n2].im = 0.0;

    for (i = 1; i < n2; ++i) {
        oup[n - i].re = oup[i].re;
        oup[n - i].im = -oup[i].im;
    }
}

auto fftC2rExec(FftRealSet* obj, FftData* inp, fft_type* oup) -> void
{
    int i;
    int n2;
    fft_type temp1;
    fft_type temp2;
    n2 = obj->cobj->N;

    auto cinp = std::make_unique<FftData[]>(n2);
    auto coup = std::make_unique<FftData[]>(n2);

    for (i = 0; i < n2; ++i) {
        temp1 = -inp[i].im - inp[n2 - i].im;
        temp2 = -inp[n2 - i].re + inp[i].re;
        cinp[i].re = inp[i].re + inp[n2 - i].re + (temp1 * obj->data[i].re) - (temp2 * obj->data[i].im);
        cinp[i].im = inp[i].im - inp[n2 - i].im + (temp2 * obj->data[i].re) + (temp1 * obj->data[i].im);
    }

    fftExec(*(obj->cobj), cinp.get(), coup.get());
    for (i = 0; i < n2; ++i) {
        oup[2 * i] = coup[i].re;
        oup[2 * i + 1] = coup[i].im;
    }
}
